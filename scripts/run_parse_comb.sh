#!/bin/bash
#SBATCH --job-name=parse_comb
#SBATCH --output=logs/parse_comb_%j.out
#SBATCH --error=logs/parse_comb_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=standard
#SBATCH --account=sbsandme_lab
#SBATCH --time=48:00:00

# split-pipe is a single-process, multi-threaded application — it uses 
# --nthreads to control parallelism internally, not MPI. So --ntasks=1 with 
# --cpus-per-task=32 is exactly right, and 
# --nthreads 32 in the split-pipe call matches the allocation.

set -euo pipefail

PROJ=/dfs9/ucightf-lab/projects/GaneA/260226_GaneA_Parse-Mega

GENOME_DIR=/dfs9/ucightf-lab/kstachel/TOOLS/REFERENCES/mouse/split-pipe/grcm38
OUTPUT_DIR=${PROJ}/output/parse_comb
NTHREADS=32

# Sublibrary 1 (already processed via snakemake parse_all rule)
# SUB1=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary1_1_S1_REP_CLEAN
SUB1=${PROJ}/output_0126I_27_GaneA_Sublibrary1_1_S1_REP_CLEAN

# Sublibraries 2-7 (from tmp; REP_CLEAN outputs from trailmaker)
SUB2=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary2_S1_REP_CLEAN
SUB3=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary3_S2_REP_CLEAN
SUB4=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary4_S1_REP_CLEAN
SUB5=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary5_S2_REP_CLEAN
SUB6=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary6_S1_REP_CLEAN
SUB7=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary7_S2_REP_CLEAN
# SUB8=${PROJ}/tmp/output_0126I_27_GaneA_Sublibrary8_REP_CLEAN  # add when available

source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate spipe

mkdir -p logs
rm -rf "${OUTPUT_DIR}"

split-pipe \
    --mode comb \
    --kit WT_mega \
    --chemistry v3 \
    --genome_dir "${GENOME_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --sublibraries "${SUB1}" "${SUB2}" "${SUB3}" "${SUB4}" "${SUB5}" "${SUB6}" "${SUB7}" \
    --nthreads "${NTHREADS}"
