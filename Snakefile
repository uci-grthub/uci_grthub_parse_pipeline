# Multi-sample CellRanger Snakemake workflow
# For processing multiple single-cell RNA-seq samples

EMAIL = "kstachel@uci.edu"
onstart:
   shell("mail -s 'STARTED' {EMAIL} < {log}")
onsuccess:
   shell("mail -s 'DONE' {EMAIL} < {log}")
onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

from datetime import datetime
import loompy
import os
import glob
import re


# Load configuration
# Load primary config and optionally a local override (config.local.yaml) if present.
import os
configfile: "config.yaml"
_local_config = "config.local.yaml"
if os.path.exists(_local_config):
    configfile: _local_config

# Extract sample list and configuration
FASTQ_DIR = config["paths"]["fastqs"]
TRANSCRIPTOME = config["references"]["transcriptome"]
OUTPUT_DIR = config["paths"]["output"]
# Optional project directory name used when building the Seurat5Shiny bundle
PROJECT_DIR_NAME = config.get("project_dir_name", "parse_comb")
REPORT = f"{OUTPUT_DIR}/{PROJECT_DIR_NAME}_report.pdf"

def create_sample_list_from_excel(excel_path, output_path):
    """
    Parse a Parse Biosciences Sample Loading Table (.xlsm) to create a
    sample_list.txt mapping sample names to their 96-well plate positions.

    Reads the 'Plate Configuration' sheet, extracts the plate layout and
    sample name table, then groups consecutive wells belonging to the same
    sample into ranges.

    Output format (one line per contiguous well block):
        {sample_name} {first_well}-{last_well}

    Sample names are lower-cased with ", " (or ",") replaced by "_".
    Single-well blocks are written as "{well}-{well}".
    """
    import openpyxl

    wb = openpyxl.load_workbook(excel_path, keep_vba=True, data_only=True)
    ws = wb["Plate Configuration"]

    PLATE_ROWS = list("ABCDEFGH")

    plate_grid = {}   # plate_row_letter -> [sample_num, ...] (12 values)
    sample_names = {} # sample_num (int) -> formatted sample name (str)

    for row in ws.iter_rows(values_only=True):
        # Plate layout rows: second cell is a letter A-H, followed by 12 ints
        if len(row) > 13 and row[1] in PLATE_ROWS:
            plate_grid[row[1]] = list(row[2:14])

        # Sample name rows: third cell is an int, fifth cell is a non-empty string
        if (len(row) > 4
                and isinstance(row[2], int)
                and isinstance(row[4], str)
                and row[4].strip()):
            sample_num = row[2]
            # Normalise: lowercase; replace spaces, ", "/","  and invalid split-pipe
            # characters "(", ")", "+" with "_", then collapse runs of "_"
            # and strip leading/trailing underscores
            sample_name = re.sub(
                r"_+", "_",
                row[4].lower()
                .replace(", ", "_")
                .replace(",", "_")
                .replace(" ", "_")
                .replace("(", "_")
                .replace(")", "_")
                .replace("+", "_")
            ).strip("_")
            sample_names[sample_num] = sample_name

    # Build ordered list of (well_name, sample_num) across all 96 wells
    wells = []
    for plate_row in PLATE_ROWS:
        for col_idx, sample_num in enumerate(plate_grid[plate_row], start=1):
            wells.append((f"{plate_row}{col_idx}", sample_num))

    # Group consecutive wells with the same sample number into runs
    lines = []
    i = 0
    while i < len(wells):
        _, sample_num = wells[i]
        j = i + 1
        while j < len(wells) and wells[j][1] == sample_num:
            j += 1
        first_well = wells[i][0]
        last_well = wells[j - 1][0]
        if sample_num is not None:
            lines.append(f"{sample_names[sample_num]} {first_well}-{last_well}")
        i = j

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return lines


# Function to retrieve sample names from config sample_list file
def get_samples_from_config(config):
    """
    Retrieve sample names from the sample_list file specified in config.
    Reads the sample_list file and extracts the first column as sample names.
    
    Args:
        config: Snakemake config dictionary
        
    Returns:
        list: Sorted list of sample names
    """
    if "sample_list" not in config:
        raise ValueError("sample_list not found in configuration")
    
    sample_list_file = config["sample_list"]
    
    if not os.path.exists(sample_list_file):
        raise FileNotFoundError(f"Sample list file not found: {sample_list_file}")
    
    samples = []
    with open(sample_list_file, 'r') as f:
        for lineno, raw_line in enumerate(f, start=1):
            # Preserve original line for error messages, but operate on a stripped copy
            line = raw_line.strip()
            if not line or line.startswith('#'):
                # Skip empty lines and comments
                continue

            # Enforce at most one space character per non-empty, non-comment line.
            # This avoids ambiguous formatting in the sample list and ensures the
            # first column is unambiguous when we split on whitespace.
            space_count = raw_line.count(' ')
            if space_count > 1:
                raise ValueError(
                    f"Invalid format in {sample_list_file} at line {lineno}:"
                    f" contains {space_count} spaces. Each non-empty line may contain at most one space."
                    f" Offending line: '{raw_line.rstrip()}'")

            # Extract first column (sample name)
            sample_name = line.split()[0] if line.split() else None
            if sample_name:
                samples.append(sample_name)
    
    if not samples:
        raise ValueError(f"No samples found in {sample_list_file}")
    
    return sorted(samples)

# Function to automatically detect sublibraries from FASTQ files
def get_sublibraries_from_fastq_dir(fastq_dir):
    """
    Automatically detect sample names from FASTQ files.
    Assumes Parse naming convention: {sample}_S{N}_L{lane}_{R1,R2,I1,I2}_001.fastq.gz
    """
    sublibraries = set()
    pattern = os.path.join(fastq_dir, "*", "*_S*_L*_R1_*.fastq.gz")
    
    for fastq_file in glob.glob(pattern):
        basename = os.path.basename(fastq_file)
        # Extract sample name before _S{number}
        match = re.match(r'^(.+)_S\d+_L\d+_R1_\d+\.fastq\.gz$', basename)
        if match:
            sample_name = match.group(1)
            sublibraries.add(sample_name)
    
    if not sublibraries:
        raise ValueError(f"No FASTQ files found in {fastq_dir} matching Parse naming pattern")
    
    return sorted(list(sublibraries))

# Generate metadata/sample_list.txt from the Excel sample loading table
_excel_files = glob.glob("metadata/*.xlsm")
if not _excel_files:
    raise FileNotFoundError("No .xlsm file found in metadata/")
create_sample_list_from_excel(_excel_files[0], "metadata/sample_list.txt")

# Automatically detect sublibraries from FASTQ directory
SUBLIBRARIES = get_sublibraries_from_fastq_dir(FASTQ_DIR)

SAMPLES = get_samples_from_config(config)

print(f"Detected sublibraries: {SUBLIBRARIES}", file=sys.stderr)
print(f"Detected samples: {SAMPLES}", file=sys.stderr)

# Optional: bp/nonbp Harmony integration (included only when the file is present)
_bp_nonbp_smk = "rules/harmony_bp_nonbp.smk"

# Rule all - defines final outputs for all samples
_bp_nonbp_outputs = (
    [
        f"{OUTPUT_DIR}/scanpy/combined_bp_harmony_integrated.h5ad",
        f"{OUTPUT_DIR}/scanpy/combined_nonbp_harmony_integrated.h5ad",
        f"{OUTPUT_DIR}/scanpy/combined_bp_harmony_integrated.rds",
        f"{OUTPUT_DIR}/scanpy/combined_nonbp_harmony_integrated.rds",
        f"{OUTPUT_DIR}/scanpy/combined_bp_raw_counts.csv",
        f"{OUTPUT_DIR}/scanpy/combined_bp_obs.csv",
        f"{OUTPUT_DIR}/scanpy/combined_bp_var.csv",
        f"{OUTPUT_DIR}/scanpy/combined_nonbp_raw_counts.csv",
        f"{OUTPUT_DIR}/scanpy/combined_nonbp_obs.csv",
        f"{OUTPUT_DIR}/scanpy/combined_nonbp_var.csv",
    ]
    if os.path.exists(_bp_nonbp_smk)
    else []
)

rule all:
    input:
        # # FastQC reports
        # expand(f"fastqc/{{sublibrary}}_r1_fastqc.html", sublibrary=SUBLIBRARIES),
        # expand(f"fastqc/{{sublibrary}}_r2_fastqc.html", sublibrary=SUBLIBRARIES),
        # CellRanger outputs
        # expand(f"{OUTPUT_DIR}/cellranger/{{sublibrary}}/outs/web_summary.html", sublibrary=SUBLIBRARIES),
        # expand(f"{OUTPUT_DIR}/cellranger/{{sublibrary}}/outs/filtered_feature_bc_matrix.h5", sublibrary=SUBLIBRARIES),
        expand(f"{OUTPUT_DIR}/parse/{{sublibrary}}/all-sample_analysis_summary.html", sublibrary=SUBLIBRARIES),
        f"{OUTPUT_DIR}/parse_comb/all_summaries.zip",
        # Parse scVI integration outputs
        # Parse Harmony integration outputs,
        # "data/merged_mm10_reference",
        # f"{OUTPUT_DIR}/seurat/parse_comb_harmony_integrated.rds",
        # f"{OUTPUT_DIR}/seurat/parse_comb_harmony_embeddings.csv",
        # f"{OUTPUT_DIR}/seurat/parse_comb_harmony_plots.pdf",
        f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb",
        f"{OUTPUT_DIR}/Seurat5Shiny/{PROJECT_DIR_NAME}",
        REPORT,
        *_bp_nonbp_outputs

# Helper function to get FASTQ files for a sample
def get_fastq_files(wildcards, read):
    """Get R1 or R2 FASTQ file for a sample"""
    pattern = f"{FASTQ_DIR}/**/{wildcards.sublibrary}_S*_L*_{read}_*.fastq.gz"
    files = glob.glob(pattern, recursive=True)
    if not files:
        raise ValueError(f"No {read} FASTQ files found for sample {wildcards.sublibrary}")
    return files

# Rule: FastQC on raw FASTQ files
rule fastqc:
    input:
        r1 = lambda wildcards: get_fastq_files(wildcards, "R1"),
        r2 = lambda wildcards: get_fastq_files(wildcards, "R2")
    output:
        r1_html = f"fastqc/{{sublibrary}}_r1_fastqc.html",
        r1_zip = f"fastqc/{{sublibrary}}_r1_fastqc.zip",
        r2_html = f"fastqc/{{sublibrary}}_r2_fastqc.html",
        r2_zip = f"fastqc/{{sublibrary}}_r2_fastqc.zip"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load fastqc/0.11.9
        mkdir -p fastqc
        fastqc -o fastqc -t {threads} {input.r1} {input.r2}
        module unload fastqc/0.11.9
        """

# Rule: MultiQC report
rule multiqc:
    input:
        expand(f"fastqc/{{sublibrary}}_r1_fastqc.html", sublibrary=SUBLIBRARIES),
        expand(f"fastqc/{{sublibrary}}_r2_fastqc.html", sublibrary=SUBLIBRARIES)
    output:
        report = f"{OUTPUT_DIR}/multiqc_report.html"
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif multiqc {OUTPUT_DIR}/fastqc -o {OUTPUT_DIR}
        module unload singularity/3.11.3
        """

# Rule: parse split-pipe for each sample
rule parse_all:
    input:
        r1 = lambda wildcards: get_fastq_files(wildcards, "R1"),
        r2 = lambda wildcards: get_fastq_files(wildcards, "R2"),
        sample_list = config["sample_list"]
    output:
        all_summaries = f"{OUTPUT_DIR}/parse/{{sublibrary}}/all-sample_analysis_summary.html",
        output_dir = directory(f"{OUTPUT_DIR}/parse/{{sublibrary}}")
    conda: "spipe"
    params:
        fastq_dir = FASTQ_DIR,
        localcores = config["params"]["localcores"],
        localmem = config["params"]["localmem"],
        parse_container = config["parse_container"],
        parse_transcriptome = config["parse_transcriptome"],
        kit = config["params"].get("kit", "WT_mega"),
        chemistry = config["params"].get("chemistry", "v3")
    threads: 8
    resources:
        mem_mb = 48000,  # 48GB in MB
        mem_mb_per_cpu=6000,
        cpus = 8,
        partition = "standard",
        account = "sbsandme_lab",
    shell:
        """
        rm -rf {output.output_dir}
        split-pipe \
        --mode all \
        --kit {params.kit} \
        --chemistry {params.chemistry} \
        --genome_dir {params.parse_transcriptome} \
        --fq1 {input.r1} \
        --fq2 {input.r2} \
        --output_dir {output.output_dir} \
        --nthreads {threads} \
        --samp_list {input.sample_list}
        """

# Rule: parse split-pipe for each sample
rule parse_comb:
    input:
        sublibraries = expand(f"{OUTPUT_DIR}/parse/{{sublibrary}}", sublibrary=SUBLIBRARIES),
        # output_dir = directory(f"{OUTPUT_DIR}/parse/{{sublibrary}}")
    output:
        all_summaries = f"{OUTPUT_DIR}/parse_comb/all_summaries.zip",
        output_dir = directory(f"{OUTPUT_DIR}/parse_comb")
    conda: "spipe"
    params:
        fastq_dir = FASTQ_DIR,
        localcores = config["params"]["localcores"],
        localmem = config["params"]["localmem"],
        parse_transcriptome = config["parse_transcriptome"],
        kit = config["params"].get("kit", "WT_mega"),
        chemistry = config["params"].get("chemistry", "v3")
    threads: 8
    resources:
        mem_mb = 131072,  # 128GB in MB
        mem_mb_per_cpu=16384,
        cpus = 8,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        rm -rf {output.output_dir}
        split-pipe \
        --mode comb \
        --kit {params.kit} \
        --chemistry {params.chemistry} \
        --genome_dir {params.parse_transcriptome} \
        --output_dir {output.output_dir} \
        --sublibraries {input.sublibraries} \
        --nthreads {threads}
        """

# Rule: create a Parse reference (mkref)
rule parse_mkref:
    input:
        fasta = config.get("mkref", {}).get("fasta", "data/merged_mm10_reference.fa"),
        genes = config.get("mkref", {}).get("genes", "data/merged_mm10_reference.gtf")
    output:
        ref_dir = directory(config.get("mkref", {}).get("output_dir", "data/merged_mm10_reference"))
    conda: "spipe"
    params:
        genome_name = config.get("mkref", {}).get("genome_name", "merged_mm10_reference")
    threads: 16
    resources:
        mem_mb = 96000,
        cpus = 16,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {output.ref_dir}
        split-pipe \
        --mode mkref \
        --genome_name {params.genome_name} \
        --fasta {input.fasta} \
        --genes {input.genes} \
        --nthreads {threads} \
        --output_dir {output.ref_dir}
        """

# Rule: Parse scVI integration
rule parse_scvi_integration:
    input:
        parse_comb_dir = f"{OUTPUT_DIR}/parse_comb"
    output:
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_scvi_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/parse_scvi_integration.py",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 4
    resources:
        mem_mb = 32000,  # 32GB in MB
        cpus = 4,
        partition = "gpu",
        account = "sbsandme_lab_gpu"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        # cd {OUTPUT_DIR}/scanpy
        python {params.script} \
            --input_dir {input.parse_comb_dir} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key}
        """

# Rule: Parse harmony integration with Seurat
rule parse_harmony_integration_r:
    input:
        parse_comb_dir = f"{OUTPUT_DIR}/parse_comb"
    output:
        seurat_obj = f"{OUTPUT_DIR}/seurat/parse_comb_harmony_integrated.rds",
        embeddings = f"{OUTPUT_DIR}/seurat/parse_comb_harmony_embeddings.csv",
        plots = f"{OUTPUT_DIR}/seurat/parse_comb_harmony_plots.pdf"
    params:
        script = "src/parse_harmony_integration.R",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = "parse_comb",
        harmony_theta = config.get("harmony_theta", 2),
        harmony_dims = config.get("harmony_dims", 30),
        cluster_resolution = config.get("cluster_resolution", 0.5),
        futures_max_size = config.get("futures_max_size", 24000)
    threads: 8
    resources:
        mem_mb = 48000,  # 48GB in MB
        cpus = 8,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        module load R/4.4.2
        mkdir -p {OUTPUT_DIR}/seurat
        Rscript {params.script} \
            --input_dir "{input.parse_comb_dir}" \
            --output_prefix "{params.output_prefix}" \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key "{params.batch_key}" \
            --output_dir "{OUTPUT_DIR}/seurat" \
            --ncores {threads} \
            --harmony_theta {params.harmony_theta} \
            --harmony_dims {params.harmony_dims} \
            --cluster_resolution {params.cluster_resolution} \
            --globals_max_size {params.futures_max_size}
        module unload R/4.4.2
        """

# Rule: Parse scVI integration
rule parse_harmony_integration_python:
    input:
        parse_comb_dir = f"{OUTPUT_DIR}/parse_comb"
    output:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined.h5ad",
        integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    conda: "scvi-tools"
    params:
        script = "src/parse_harmony_integration.py",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined"
    threads: 16
    resources:
        mem_mb = 240000,  # 32GB in MB
        cpus = 16,
        partition = "hugemem",
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        # cd {OUTPUT_DIR}/scanpy
        python {params.script} \
            --input_dir {input.parse_comb_dir} \
            --output_prefix {params.output_prefix} \
            --min_genes {params.min_genes} \
            --min_cells {params.min_cells} \
            --n_top_genes {params.n_top_genes} \
            --batch_key {params.batch_key}
        """

# Rule: parse harmony notebook
rule parse_harmony_notebook:
    input:
       integration_results = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    output:
        f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb"
    conda: "scvi-tools"
    params:
        script = "src/submit_harmony_integration.sh",
        input_dir = f"{OUTPUT_DIR}/cellranger",
        min_genes = config.get("min_genes", 300),
        min_cells = config.get("min_cells", 5),
        n_top_genes = config.get("n_top_genes", 2000),
        batch_key = config.get("batch_key", "batch"),
        output_prefix = f"{OUTPUT_DIR}/scanpy/combined",
        groupby_var = config.get("condition", "condition"),
        group1_value = config.get("group1_value", "treat"),
        group2_value = config.get("group2_value", "control"),
        notebook = config["notebooks"]["harmony_integration_notebook"]
    threads: 4
    resources:
        mem_mb = 32000,  # 32GB in MB
        cpus = 8,
        account = "sbsandme_lab"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scanpy
        echo {input.integration_results}
        {params.script} --no-integration \
        --min-genes {params.min_genes} \
        --groupby_var {params.groupby_var} \
        --group1_value {params.group1_value} \
        --group2_value {params.group2_value} \
        --notebook {params.notebook} \
        {params.output_prefix}
        """

rule save_raw_adata:
    input:
        combined_adata = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    output:
        adata = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated_data.h5ad",
        scale_data_adata = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated_scale_data.h5ad"
    threads: 20
    params:
        script = "src/save_raw_adata.py"
    resources:
        mem_mb = 320000,  # 300GB in MB
        partition = "hugemem",
        cpus = 20,
        account = "sbsandme_lab"
    shell:
        '''
        {params.script} {input.combined_adata} {output.adata} {output.scale_data_adata}
        '''

rule convert_harmony_integrated_h5ad_to_rds:
    input:
        adata = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.h5ad"
    output:
        rds = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.rds"
    params:
        script = "src/convert_h5ad_to_seurat.R",
        output_prefix = f"scaled_data"
    threads: 20
    resources:
        mem_mb = 320000,  # 300GB in MB
        partition = "hugemem",
        cpus = 20,
        account = "sbsandme_lab"
    shell:
        '''
        module load R/4.3.3 
        Rscript {params.script} --data {input.adata} --output {output.rds}
        module unload R/4.3.3
        '''

rule build_seurat5shiny:
    input:
        rds = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated.rds",
        scale_data_adata = f"{OUTPUT_DIR}/scanpy/combined_harmony_integrated_scale_data.h5ad"
    output:
        shiny_dir = directory(f"{OUTPUT_DIR}/Seurat5Shiny/{PROJECT_DIR_NAME}"),
        rds = f"{OUTPUT_DIR}/Seurat5Shiny/{PROJECT_DIR_NAME}/seurat5.rds"
    params:
        app_title = config.get("shiny_app", {}).get("title", PROJECT_DIR_NAME),
        r_script = "src/seurat_scale_data.R",
    threads: 20
    resources:
        mem_mb = 320000,  # 300GB in MB
        partition = "hugemem",
        cpus = 20,
        account = "sbsandme_lab"
    shell:
        '''
        module load R/4.3.3 
        mkdir -p {OUTPUT_DIR}/Seurat5Shiny
        cp -r /dfs9/ucightf-lab/kstachel/Seurat5Shiny {output.shiny_dir}
        echo {params.app_title} > {output.shiny_dir}/title.txt
        date > {output.shiny_dir}/restart.txt
        {params.r_script} {input.rds} {output.rds}  
        module unload R/4.3.3 
        '''

# Rule: Generate multi-sample summary
rule multi_sample_summary:
    input:
        web_summaries = expand(f"{OUTPUT_DIR}/parse/{{sublibrary}}/outs/web_summary.html", sublibrary=SUBLIBRARIES)
    output:
        summary = f"{OUTPUT_DIR}/multi_sample_summary.txt"
    run:
        with open(output.summary, 'w') as f:
            f.write("CellRanger Multi-Sample Analysis Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total samples processed: {len(SUBLIBRARIES)}\n\n")
            
            for sample in SUBLIBRARIES:
                f.write(f"Sample: {sample}\n")
                f.write(f"  - Web summary: {OUTPUT_DIR}/cellranger/{sample}/outs/web_summary.html\n")
                f.write(f"  - Filtered matrix: {OUTPUT_DIR}/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5\n")
                f.write(f"  - BAM file: {OUTPUT_DIR}/cellranger/{sample}/outs/possorted_genome_bam.bam\n\n")

# Rule: loompy
rule loompy:
    input:
        folder = f"{OUTPUT_DIR}/cellranger/{{sublibrary}}"
    output:
        loom_file = f"{OUTPUT_DIR}/loom/{{sublibrary}}.loom"
    resources:
        mem_mb = 8000,  # 8GB for loompy
        cpus = 4,
        partition = "standard"
    threads: 4
    run:
        loompy.create_from_cellranger(input.folder, f"{OUTPUT_DIR}/loom", genome="hg38")

# Rule: Seurat
rule seurat:
    input:
        matrix_dir = f"{OUTPUT_DIR}/{{sublibrary}}/outs/filtered_feature_bc_matrix/",
        script = "src/process_seurat.R"
    output:
        seu_path = f"{OUTPUT_DIR}/seurat/{{sublibrary}}_seu.rds"
    log:
        f"{OUTPUT_DIR}/logs/seurat_{{sublibrary}}.log"
    threads: 4
    shell:
        '''Rscript {input.script} matrix_dir="{input.matrix_dir}" seu_path="{output.seu_path}" celltype_ref="config['celltype_ref']" > {log} 2>&1'''

# Rule: Save Seurat embeddings
rule seurat_embeddings:
    input:
        seu_path = f"{OUTPUT_DIR}/seurat/{{sublibrary}}_seu.rds",
        nb_path = f"{OUTPUT_DIR}/numbat/{{sublibrary}}_numbat.rds",
        script = "src/save_seurat_embeddings.R"
    output:
        seurat_embeddings = f"{OUTPUT_DIR}/seurat/{{sublibrary}}_embeddings.csv"
    log:
        f"{OUTPUT_DIR}/logs/seurat_{{sublibrary}}_embeddings.log"
    threads: 2
    shell:
        '''Rscript {input.script} seu_path="{input.seu_path}" nb_path="{input.nb_path}" > {log} 2>&1'''

# Rule: Velocyto
rule velocyto:
    input:
        sample_folder = f"{OUTPUT_DIR}/{{sublibrary}}"
    output:
        loom = f"{OUTPUT_DIR}/velocyto/{{sublibrary}}.loom"
    threads: 4
    shell:
        """
        module load velocyto/0.17.17
        velocyto run10x -m config["references"]["repeat_mask"] {input.sample_folder} config["gtf"] > {output.loom}
        module unload velocyto/0.17.17
        """

# Rule: scVelo
rule scvelo:
    input:
        loom_file = f"{OUTPUT_DIR}/velocyto/{{sublibrary}}.loom",
        anndata_file = f"{OUTPUT_DIR}/scanpy/{{sublibrary}}.h5ad",
        seurat_embeddings = f"{OUTPUT_DIR}/seurat/{{sublibrary}}_embeddings.csv"
    output:
        scvelo_h5ad = f"{OUTPUT_DIR}/scanpy/{{sublibrary}}_scvelo.h5ad"
    threads: 2
    shell:
        '''python src/compute_velocity.py {input.anndata_file} {input.loom_file} {input.seurat_embeddings}'''

# Rule: pileup and phasing
rule pileup_and_phasing:
    input:
        bam = f"{OUTPUT_DIR}/{{sublibrary}}/outs/possorted_genome_bam.bam",
        barcodes = f"{OUTPUT_DIR}/{{sublibrary}}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        script = "src/pileup_and_phase.R"
    output:
        allele_df = f"{OUTPUT_DIR}/numbat/{{sublibrary}}_allele_counts.tsv.gz"
    threads: 4
    shell:
        '''Rscript {input.script} --label {wildcards.sublibrary} --samples {wildcards.sublibrary} --bams {input.bam} --barcodes {input.barcodes} --gmap config["gmap"] --snpvcf config["snpvcf"] --paneldir config["paneldir"] --outdir {OUTPUT_DIR}/numbat/{wildcards.sublibrary} --ncores {threads}'''

# Rule: Numbat
rule numbat:
    input:
        allele_df = f"{OUTPUT_DIR}/numbat/{{sublibrary}}_allele_counts.tsv.gz",
        matrix_file = f"{OUTPUT_DIR}/{{sublibrary}}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        seu_path = f"{OUTPUT_DIR}/seurat/{{sublibrary}}_seu.rds",
        script = "src/run_numbat.R"
    output:
        done_file = f"{OUTPUT_DIR}/numbat/{{sublibrary}}/done.txt"
    threads: 4
    shell:
        '''Rscript {input.script} seu_path="{input.seu_path}" ref_path=config["ref_path"] tau=config["tau"] read_prop=config["read_prop"] max_iter=config["max_iter"] min_LLR=config["min_LLR"] t=config["numbat_t"] cell_ceiling=config["cell_ceiling"] max_entropy=config["max_entropy"] allele_df="{input.allele_df}" matrix_file="{input.matrix_file}" out_dir="{OUTPUT_DIR}/numbat/{wildcards.sublibrary}" ncores={threads} rprof_out={OUTPUT_DIR}/numbat/{wildcards.sublibrary}/log.prof > {output.done_file} 2>&1'''

# Rule: SCENIC analysis
include: "rules/scenic.smk"

# Rule: CellBender remove-background
rule cellbender:
    input:
        h5 = f"{OUTPUT_DIR}/cellranger/{{sublibrary}}/outs/raw_feature_bc_matrix.h5"
    output:
        clean_h5 = f"{OUTPUT_DIR}/cellbender/{{sublibrary}}_filtered.h5"
    conda:
        "envs/cellbender_env.yaml"
    threads: 8
    resources:
        mem_mb = 64000,
        cpus = 8,
        partition = "standard"
    params:
        expected_cells = config.get("cellbender_expected_cells", 3000),
        total_droplets = config.get("cellbender_total_droplets", 20000)
    shell:
        """
        cellbender remove-background \
            --input {input.h5} \
            --output {output.clean_h5} \
            --expected-cells {params.expected_cells} \
            --total-droplets-included {params.total_droplets} \
            --cpu-threads {threads}
        """

# Rule: Scrublet doublet detection
rule scrublet:
    input:
        matrix = f"{OUTPUT_DIR}/cellranger/{{sublibrary}}/outs/filtered_feature_bc_matrix.h5"
    output:
        doublets = f"{OUTPUT_DIR}/scrublet/{{sublibrary}}_doublets.csv"
    threads: 2
    resources:
        mem_mb = 8000,
        cpus = 2,
        partition = "standard"
    params:
        script = "src/run_scrublet.py"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/scrublet
        python {params.script} --input {input.matrix} --output {output.doublets}
        """

rule generate_report:
    input:
        parse_comb = f"{OUTPUT_DIR}/parse_comb/all_summaries.zip",
        notebook = f"{OUTPUT_DIR}/scanpy/inspect_integrated_anndata_combined.ipynb",
        metadata = "metadata/metadata.csv",
        script = "src/generate_report.py"
    output:
        report = REPORT
    conda: "scvi-tools"
    params:
        parse_dir = f"{OUTPUT_DIR}/parse_comb",
        fastq_dir = FASTQ_DIR
    threads: 1
    resources:
        mem_mb = 4000,
        cpus = 1,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        pip install --quiet reportlab 2>/dev/null || true
        python {input.script} \
            --fastq-dir {params.fastq_dir} \
            --parse-dir {params.parse_dir} \
            --metadata {input.metadata} \
            --output {output.report}
        """

if os.path.exists(_bp_nonbp_smk):
    include: _bp_nonbp_smk
