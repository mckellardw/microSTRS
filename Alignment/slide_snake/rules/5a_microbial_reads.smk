#####################################
# After running standard cell ranger RNA, extract unmapped reads from possorted_bam_unmapped
#####################################
# extract unmapped reads and run Kraken on unmapped reads
rule getUnmappedReads:
	input:
		'{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
	params:
		krakenDB=config['KRAKEN_DB']
	output:
		one='{OUTDIR}/{sample}/STARsolo/possorted_bam_kraken2.out.gz',
		two='{OUTDIR}/{sample}/STARsolo/kraken'
	threads: CORES
	shell:
		"""
		{KRAKEN} --use-names --db {params.krakenDB} --minimum-hit-groups 3 --report {output.two} --use-mpa-style <(samtools view -b -f 4 {input} | samtools fasta) | gzip > {output.one}
		"""

# extract cell barcode and UMI of unmapped reads
rule unmappedBarcodes:
	input:
		bam='{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam',
		kraken='{OUTDIR}/{sample}/STARsolo/possorted_bam_kraken2.out.gz'
	output:
		'{OUTDIR}/{sample}/STARsolo/Kraken_barcodes.txt.gz'
	threads: CORES
	shell:
		"""
		paste <(zcat {input.kraken}) <(samtools view -f 4 {input.bam} | grep -o -P '(?<=CB:Z:).*(?=UB:Z:)') <(samtools view -f 4 {input.bam} | grep -o -P '(?<=UB:Z:).*') | gzip > {output}
		"""
rule cache_preQC_h5ad_STAR_raw:
    input:
        BCS = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz',
        GENES = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz',
        MAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx.gz',
        BB_map = lambda wildcards: BB_DICT[wildcards.sample]
    output:
        H5AD = '{OUTDIR}/{sample}/STARsolo/Microbial.h5ad'
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx()
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_mtx_to_h5ad.py \
            --mat_in {input.MAT} \
            --feat_in {input.GENES} \
            --bc_in {input.BCS} \
            --bb_map {input.BB_map}\
            --ad_out {output.H5AD}\
            --feat_col 1 \
            --remove_zero_features
            """
        )
