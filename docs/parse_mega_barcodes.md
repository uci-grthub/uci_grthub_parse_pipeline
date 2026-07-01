The 4th round (sub-library PCR) uses the Illumina sequencing index — it doesn't appear as a biological barcode in R2. Split-pipe handles that by processing each sub-library's FASTQs separately (--mode all) then combining (--mode comb).

The n302 barcode count now makes sense: each of the 96 wells has two barcode sequences:

stype	Count	Meaning
R	96	Random hexamer RT primer barcode (one per well)
T	96	Poly-dT RT primer barcode (one per well)
X_*	~110	Immune profiling special barcodes
Parse uses both oligo-dT and random hexamer primers simultaneously in the RT step. A cell's transcript gets whichever primer type captured it, and split-pipe records this as the ptype field in the preprocessed FASTQ read header (the __T__ or __R__ field in read names like 20_15_07__T__20_63_55__...).

So the corrected summary:

3 biological barcode rounds are captured in R2 (not 4)
Round 3 / bc3 (R3_v3, 96 barcodes): RT barcode — the sample-defining step
Round 2 / bc2 (v1, 96 barcodes): first ligation
Round 1 / bc1 (n302_R1_v3_8, 302 entries): second ligation — 96 wells × 2 primer types (R + T) = 192 barcodes plus ~110 immune extras
4th "round" (sub-library PCR): the Illumina sequencing index — handled upstream by bcl2fastq/demultiplexing, not by split-pipe barcode matching
The whitelist and linker sequences I listed earlier are still correct; the complication is just that bc1 has duplicate sequences-per-well (R vs T) rather than 302 independent wells.
