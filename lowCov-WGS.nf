#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// INPUT CHANNELS
    // Metadata and short-read sequencing reads (paired-end or single-end) as specified 
    // in a CSV
    ch_reads = Channel
        .fromPath( params.samplesheet )
        .splitCsv( header: true )
        .map { row -> tuple(row.sample, row.population, row.species, row.library_prep, row.seq_platform, file(row.reads1_path), file(row.reads2_path)) }
	
	ch_reference = Channel
		.fromPath( params.reference )
	
	// WORKFLOW STEPS
    // Routine bioinformatic processing and QC
    INTERLEAVE_READS (
		ch_reads
	)

	REMOVE_OPTICAL_DUPLICATES (
		INTERLEAVE_READS.out
	)

	REMOVE_LOW_QUALITY_REGIONS (
		REMOVE_OPTICAL_DUPLICATES.out
	)

	TRIM_ADAPTERS (
		REMOVE_LOW_QUALITY_REGIONS.out
	)

	REMOVE_ARTIFACTS (
		TRIM_ADAPTERS.out
	)

	ERROR_CORRECT_PHASE_ONE (
		REMOVE_ARTIFACTS.out
	)

	ERROR_CORRECT_PHASE_TWO (
		ERROR_CORRECT_PHASE_ONE.out
	)

	ERROR_CORRECT_PHASE_THREE (
		ERROR_CORRECT_PHASE_TWO.out
	)

	NORMALIZE_READS (
		ERROR_CORRECT_PHASE_THREE.out
	)

	MERGE_READS (
		NORMALIZE_READS.out
	)

	QUALITY_TRIM (
		MERGE_READS.out
	)

    // ORIENT_READS (
    //     QUALITY_TRIM.out,
	// 	ch_reference
    // )

    FASTP_FILTER (
        QUALITY_TRIM.out
    )

    // FASTQC (
    //     FASTP_FILTER.out
    // )

    // MULTIQC (
    //     FASTQC.out.collect()
    // )

    MAP_TO_REFERENCE (
        FASTP_FILTER.out,
		ch_reference
    )

    // ASSESS_DEPTH (
    //     MAP_TO_REFERENCE.out
    // )

    // PLOT_DEPTH (
    //     ASSESS_DEPTH.out
    // )

    // ANGSD_GL (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> sample, population, species, library_prep, bam }
    //         .groupTuple( by:[2,3] )
    // )

    // CALL_VARIANTS (
    //     MAP_TO_REFERENCE.out
    // )

    // FILTER_VARIANTS (
    //     CALL_VARIANTS.out
    // )

    // MERGE_VARIANTS (
    //     FILTER_VARIANTS.out.vcf
	// 		.map { species, library_prep, vcf -> tuple( file(raw_vcf), species, prep_type ) }
	// 		.filter { file(it[0]).countLines() > 0 }
	// 		.groupTuple( by: [1,2] ),
	// 	FILTER_VARIANTS.out.index.collect()
    // )

	// FILTER_SNPS_BY_SAMPLE_COUNT (
	// 	MERGE_VARIANTS.out.vcf
	// )

	// RECORD_SNP_FILTERS (
	// 	FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf
	// )


    // Individual-level analyses
    // STRUCTURE (
    //     MERGE_VARIANTS.out.vcf
    // )

    // NGSADMIX (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
    // )

    // PCANGSD (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
	// )

    // FASTPCA (
	// 	ANGSD_GL.out
	// )

    // OHANA (
    //     MERGE_VARIANTS.out.vcf
	// )

    // ANGSD_GWAS (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
	// )

    // NGSLD (
	// 	FILTER_VARIANTS.out.vcf
	// )

    // ROH (
	// 	MAP_TO_REFERENCE.out
	// )

	// PEDIGREE_STRUCTURES (
	// 	MERGE_VARIANTS.out.vcf
	// )

    // D_STATISTIC ()

    // Population-level analyses
    // ANGSD_AF (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
	// )
    
    // ANGSD_SFS (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
	// )

    // VISUALIZE_SFS (
	// 	ANGSD_SFS.out.sfs
	// )

    // STAIRWAY_PLOT (
	// 	ANGSD_SFS.out.sfs
	// )

    // ANGSD_FST (
    //     MAP_TO_REFERENCE.out
    //         .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
    //         .groupTuple( by:[0,1] )
	// )

    // NGSTOOLS_FST ()

    // VCFLIB_ASSESS_FST (
	// 	ANGSD_FST.out
	// 		.mix (
	// 			NGSTOOLS_FST.out
	// 		)
	// )
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// specifying whether to run in low disk mode
if( params.low_disk_mode == true ) {
	params.publishMode = 'symlink'
}
else {
	params.publishMode = 'copy'
}

// Preprocessing results subdirectories
params.preprocessing = params.results + "/01_preprocessing"
params.optical_dedupe = params.preprocessing + "/01_optical_dedup"
params.low_quality = params.preprocessing + "/02_remove_low_quality"
params.trim_adapters = params.preprocessing + "/03_trim_adapters"
params.remove_artifacts = params.preprocessing + "/04_remove_artifacts"
params.error_correct = params.preprocessing + "/05_error_correct"
params.normalize = params.preprocessing + "/06_normalized_reads"
params.merged_reads = params.preprocessing + "/07_read_merging"
params.qtrim = params.preprocessing + "/08_quality_trim"
// params.orient = params.preprocessing + "/09_orient_reads"
params.fastp = params.preprocessing + "/09_fastp_filter"

// FASTQ read QC reports
params.read_reports = params.preprocessing + "/10_read_reports"
params.fastqc = params.read_reports + "/FastQC"
params.multiqc = params.read_reports + "/MultiQC"

// Read mapping, SNP calling, and SNP filtering
params.read_mapping = params.results + "/02_read_mapping"
params.variant_calling = params.results + "/03_snp_calling"
params.per_sample = params.variant_calling + "/by_sample_snps"
params.merged = params.variant_calling + "/merged_snps"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

// Standard bioinformatic processing and QC. 
// -----------------------------------------
// These steps are designed for Illumina short reads and not to any other 
// platform or read length, but support for these sequence read configurations 
// may be added in the future. The reads preprocessing steps were inspired by
// standards at the Joint Genomics Institute, as described at:
// https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/data-preprocessing/

process INTERLEAVE_READS {
	
	/* 
	First, we interleave paired reads into a single FASTQ, or give 
	single-end reads a more informative name. This step is merely
	a convenience, and does not alter or remove any reads in the
	raw FASTQ files.
	*/

	tag "${sample}"
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path(reads2_path)
	
	output:
    tuple val(sample), val(population), val(species), val(library_prep), env(paired_status), path("*.fastq.gz")
	
	script:
    if ( reads2_path.exists() ){
        """
		reformat.sh \
		in1=`realpath ${reads1_path}` \
		in2=`realpath ${reads2_path}` \
		out=${sample}_${species}_${library_prep}.fastq
		paired_status=`echo "paired"`
        """
    } else {
        """
        cp `realpath ${reads1_path}` ./${sample}_${species}_${library_prep}_se.fastq.gz
		paired_status=`echo "single-end"`
        """
    }
	
}


process REMOVE_OPTICAL_DUPLICATES {

	/* 
	This process removes optical duplicates from the Illumina flow cell.
	*/

	tag "${sample}"
	publishDir params.optical_dedupe, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	clumpify.sh in=`realpath ${reads}` \
	out=${sample}_clumped.fastq.gz \
	threads=${task.cpus} \
	dedupe optical tossbrokenreads
	"""

}


process REMOVE_LOW_QUALITY_REGIONS {

	/* 
	Low quality regions of each read are removed in this process.
	*/

	tag "${sample}"
	publishDir params.low_quality, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	filterbytile.sh in=`realpath ${reads}` \
	out=${sample}_filtered_by_tile.fastq.gz \
	threads=${task.cpus}
	"""

}


process TRIM_ADAPTERS {
	
	/* 
	This process takes the pipeline's first pass at removing adapters, along
	with some additional parameters recommended by bbmap.
	*/

	tag "${sample}"
	publishDir params.trim_adapters, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	if ( paired_status == "paired" )
		"""
		bbduk.sh in=`realpath ${reads}` \
		out=${sample}_trim_adapters.fastq.gz \
		ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered \
		threads=${task.cpus}
		"""
	else
		"""
		bbduk.sh in=`realpath ${reads}` \
		out=${sample}_trim_adapters.fastq.gz \
		ktrim=r k=23 mink=11 hdist=1 minlen=70 ref=adapters ftm=5 ordered \
		threads=${task.cpus}
		"""

}


process REMOVE_ARTIFACTS {

	/* 
	Here we remove various contantimants that may have ended up in the reads,
	such as PhiX sequences that are often used as a sequencing control.
	*/

	tag "${sample}"
	publishDir params.remove_artifacts, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	bbduk.sh in=`realpath ${reads}` \
	out=${sample}_remove_artifacts.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_ONE {

	/* 
	Bbmap recommends three phases of read error correction, the first of which
	goes through BBMerge.
	*/

	tag "${sample}"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	bbmerge.sh in=`realpath ${reads}` \
	out=${sample}_error_correct1.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample}_ihist_merge1.txt \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_TWO {

	/* 
	The second phase of error correction goes through clumpify.sh
	*/

	tag "${sample}"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	clumpify.sh in=`realpath ${reads}` \
	out=${sample}_error_correct2.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_THREE {

	/* 
	The third phase of error correction uses tadpole.sh.
	*/

	tag "${sample}"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	tadpole.sh in=`realpath ${reads}` \
	out=${sample}_error_correct3.fastq.gz \
	ecc k=62 ordered \
	threads=${task.cpus}
	"""

}


process NORMALIZE_READS {

	/* 
	In the process, reads are "normalized" by removing exceptionally repetitive
	sequences. Sequences that repeat above a threshold number of times (100) 
	are considered outliers and removed from sequences with "normal" numbers
	of repetitions in the dataset.
	*/

	tag "${sample}"
	publishDir params.normalize, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	bbnorm.sh in=`realpath ${reads}` \
	out=${sample}_normalized.fastq.gz \
	target=100 \
	hist=${sample}_khist.txt \
	peaks=${sample}_peaks.txt
	"""

}


process MERGE_READS {
	
	/* 
	Next, we merge any perfectly identical reads to reduce computational overhead
	without losing read support.
	*/

	tag "${sample}"
	publishDir params.merged_reads, mode: params.publishMode, overwrite: true

	errorStrategy 'ignore'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)
	
	output:
    tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*_merged.fastq.gz")
	
	script:
	if ( paired_status == "paired" )
		"""
		bbmerge-auto.sh in=`realpath ${reads}` \
		out=${sample}_merged.fastq.gz \
		outu=${sample}_unmerged.fastq.gz \
		strict k=93 extend2=80 rem ordered \
		ihist=${sample}_ihist_merge.txt \
		threads=${task.cpus}
		"""
	else
		"""
		cp `realpath ${reads}` ./${sample}_se_not_merged.fastq.gz
		"""
	
}


process QUALITY_TRIM {

	/* 
	Here we quality trim reads from both ends to a minimum Phred quality of 10, 
	and enforce a minimum read length of 70 bases. 
	*/

	tag "${sample}"
	publishDir params.qtrim, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")

	script:
	"""
	bbduk.sh in=`realpath ${reads}` \
	out=${sample}_qtrimmed.fastq.gz \
	qtrim=rl trimq=10 minlen=70 ordered \
	threads=${task.cpus}
	"""

}


// process ORIENT_READS {
	
// 	/* 
// 	Next, reads are oriented such that they have the same 5' to 3' polarity as the 
// 	reference sequence, which reduces the need to map or primer-trim reverse 
// 	complements. This step may be removed in future versions.
// 	*/
	
// 	tag "${sample}"
// 	publishDir params.orient, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

// 	errorStrategy 'ignore'
	
// 	input:
//     tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)
// 	each path(reference)
	
// 	output:
// 	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")
	
// 	script:
// 	"""
//     vsearch --orient `realpath ${reads}` \
// 	--db ${reference} \
//     --fastqout ${sample}_oriented.fastq.gz
// 	"""
// }


process FASTP_FILTER {
	
	/* 
	To be safe, we use a second algorithm here to trim adapters, quality-
	control, and generate a before and after report. Fastp also offers
	a low-complexity trimmer, which can be useful for short reads.
	*/
	
	tag "${sample}"
	publishDir params.fastp, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy 'ignore'
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.fastq.gz")
	
	script:
	if ( paired_status == "paired" )
		"""
		fastp --in1 `realpath ${reads}` \
		--out1 ${sample}_${species}_${library_prep}_filtered.fastq.gz \
		--qualified_quality_phred 15 \
		--length_required 70 \
		--detect_adapter_for_pe --correction \
		--low_complexity_filter \
		--trim_tail1 5 \
		--thread ${task.cpus} \
		--html ${sample}_${species}_${library_prep}.html
		"""
	else
		"""
		fastp --in1 `realpath ${reads}` \
		--out1 ${sample}_filtered.fastq.gz \
		--qualified_quality_phred 15 \
		--length_required 70 \
		--low_complexity_filter \
		--trim_tail1 5 \
		--thread ${task.cpus} \
		--html ${sample}_${species}_${library_prep}.html
		"""
}


process FASTQC {
	
	/* 
	Read preprocessing is now complete. To assess the quality of your
	dataset, we first run FastQC on each sample's reads.
	*/
	
	tag "${sample}"
	publishDir params.fastqc, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy 'ignore'
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)
	
	output:
	path "*"
	
	script:
	"""
	fastqc -f fastq -o ./ -t ${task.cpus} --extract ${reads}
	"""
}


process MULTIQC {
	
	/* 
	Finally, we collate the individual FastQC reports into one MultiQC report.
	*/
	
	publishDir params.multiqc, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy 'ignore'
	
	cpus 4
	
	input:
	path fastqc_files
	
	output:
	path "*"
	
	script:
	"""
	multiqc --threads ${task.cpus} .
	"""
}


process MAP_TO_REFERENCE {
	
	/* 
	This step uses the BWA MEM algorithm to map each FASTQ to a reference genome.
	It then uses samtools to convert the resulting SAM file into a BAM file, which
	is binary and thus much smaller.
	*/
	
	tag "${sample}"
	publishDir params.read_mapping, mode: 'copy', overwrite: true

	errorStrategy 'ignore'
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(reads)
	each path(reference)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path("*.bam")
	
	script:
	"""
	bwa index ${reference} && \
	bwa mem -t ${task.cpus} ${reference} `realpath ${reads}` \
    | samtools view -bS - > ${sample}_${species}_${library_prep}.bam
	"""
}


process ASSESS_DEPTH {
	
	/* 
	Next we generate depth reports for each BAM file using mosdepth.
	*/
	
	publishDir params.depth, mode: 'copy'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(paired_status), path(bam)
	
	output:
	path "*"
	
	script:
	"""
	mosdepth -t ${task.cpus} ${sample}_${species}_${library_prep} ${bam}
	"""
}


process PLOT_DEPTH {
	
	/* 
	The mosdepth results are then plotted with a simple R script.
	*/
	
	publishDir params.depth, mode: 'copy'
	
	input:
	path mosdepth_files
	
	output:
	path "*png"
	
	script:
	"""
	plot_depth_per_scaffold.R
	"""
}


process ANGSD_GL {
	
	/* 
	BAM files are then used as inputs for ANGSD, which will generate genotype
	likelihoods that can be used more confidently with low coverage data than
	hard variant calls.
	*/
	
	tag "${species} ${library_prep}"
	publishDir params.gl, mode: 'copy'
	
	cpus 8
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:


	when:
	params.variant_call_only == false
	
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	angsd -bam bam_file_list.txt \
    -ref ${params.reference} \
    -doMajorMinor 1 -doMaf 1 -doGlf 2 \
    -minInd ${params.min_samples} \
    --nThreads ${task.cpus} \
    -out ${species}_${prep_type}
	"""
}


process CALL_VARIANTS {
	
	/* 
	This process does something described here
	*/
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:
	tuple val(species), val(library_prep), path("*.vcf.gz")
	
	script:
	"""
	
	"""
}


process FILTER_VARIANTS {
	
	/* 
	This process does something described here
	*/
	
	tag "${species} ${library_prep}"
	publishDir params.per_sample, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(vcf)
	
	output:
	tuple val(species), val(library_prep), path("*.vcf.gz"), emit: vcf
	path "*.tbi", emit: index
	
	script:
	"""
	
	vcftools --gzvcf ${vcf} \
	--max-alleles 2 \
	--mac ${params.minor_allele_count} \
	--max-missing-count ${params.max_snp_missingness} \
	--minQ ${params.min_quality} \
	--minDP ${params.min_depth} \
	--remove-indels \
	--remove-filtered-all \
	--recode --stdout \
	| bgzip -c > "${sample}_${prep}_filtered.vcf.gz" && \
	tabix -p vcf "${sample}_${prep}_filtered.vcf.gz"
    
	"""
}


process MERGE_VARIANTS {
	
	/* 
	This process does something described here
	*/
	
	tag "${species} ${library_prep}"
	
	input:
	tuple path(vcf_files), val(species), val(prep)
	path index_files
	
	output:
	tuple path("*.vcf.gz"), val(species), val(prep), env(sample_size), emit: vcf
	path "*.txt"
	
	shell:
	'''
	
	bcftools merge \
	--merge snps \
	--output-type z \
	--threads !{task.cpus} \
	--output !{species}_multisample_!{prep}.vcf.gz \
	*!{species}*.vcf.gz
	
	'''
}


process FILTER_SNPS_BY_SAMPLE_COUNT {
	
	tag "${prep} ${species}"
	publishDir params.merged, pattern: "*.vcf.gz", mode: 'copy'
	publishDir params.merged, pattern: "*.tbi", mode: 'copy'
	
	cpus 2
	
	input:
	tuple path(vcf), val(species), val(prep)
	
	output:
	tuple path("${species}_${prep}_filtered.vcf.gz"), val(species), val(prep), env(sample_size), emit: vcf
	path "*.tbi", emit: index
	
	script:
	"""
	
	vcftools --gzvcf ${vcf} \
	--max-alleles 2 \
	--mac ${params.minor_allele_count} \
	--max-missing-count ${params.max_snp_missingness} \
	--remove-indels \
	--remove-filtered-all \
	--recode --stdout \
	| bgzip -c > "${species}_${prep}_filtered.vcf.gz" && \
	tabix -p vcf "${species}_${prep}_filtered.vcf.gz" && \
	sample_size=`bcftools query -l ${species}_${prep}_filtered.vcf.gz | wc -l`
	
	"""
	
}


process RECORD_SNP_FILTERS {
	
	tag "${prep} ${species}"
	publishDir params.variant_calling, mode: 'copy', pattern: '*.txt'
	
	cpus 1
	
	input:
	tuple path(vcf), val(species), val(prep), val(sample_size)
	
	output:
	path "*.txt"
	
	shell:
	'''
	sample_ids=`bcftools query -l !{vcf}`
	minInd="$((!{sample_size} - !{params.max_snp_missingness}))"
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele count: SNP must be observed in !{params.minor_allele_count} out of !{sample_size} samples" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / !{sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: " >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo ${sample_ids} >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	'''
	
}



// Individual-level analyses

process STRUCTURE {
	
	/* 
	This process does something described here
	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple path(merged_vcf), val(species), val(library_prep), val(sample_size)
	
	output:
	path "*"
	
	when:
	params.estimate_structure == true && params.variant_call_only == false
	
	script:
	"""
    plink --vcf ${merged_vcf} --recode structure \
    --out ${species}_${library_prep} && \
	structure -K 1:10 -o . \
    --input=${species}_${library_prep}.str \
    --format=str --full \
    --seed=${params.random_seed}
	"""

}


process NGSADMIX {
	
	/* 
	This process does something described here
	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(bam_files)
	
	output:
	path "*"

	when:
	params.variant_call_only == false
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	ngsAdmix -P bam_file_list.txt \
    -B ${params.reference} \
    -K 3 \
    -o ${species}_${library_prep}
	"""
}


process PCANGSD {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(species), val(library_prep), path(beagle)
	
	output:
	

	when:
	params.variant_call_only == false && params.principal_component_analyses == true
	
	script:
	"""
	PCAngsd -beagle ${beagle} \
	-o ${species}_${library_prep} \
	-minMaf ${params.minor_allele_frequency} \
	-admix \
	-threads ${task.cpus}
	"""
}


process FASTPCA {
	
	/* 

	Important note to deal with from ChatGPT:

	"Note that fastpca assumes that the genotype data is in Hardy-Weinberg 
	equilibrium (HWE). If your data violates HWE assumptions, the PCA results 
	may be biased. You may need to filter out SNPs that deviate from HWE before 
	running fastpca.

	Also note that fastpca assumes that the genotype data is diploid and biallelic. 
	If your data includes polyploid or multiallelic sites, you may need to preprocess 
	the data to convert it to a diploid, biallelic format."

	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(species), val(library_prep), path(vcf)
	
	
	output:
	
	
	when:
	params.variant_call_only == false && params.principal_component_analyses == true
	
	script:
	"""
	plink --vcf input_file.vcf \
	--make-bed --out ${species}_${library_prep} && \
	fastpca -b ${species}_${library_prep}.bed \
	-f ${species}_${library_prep}.bim \
	-o ${species}_${library_prep} \
	-t ${task.cpus} \
	-n 10
	"""
}


process OHANA {
	
	/* 

	NOTE: 

	Need to work out how to use Plink and another program to create the input file
	formats that ohana requires.

	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'

	cpus 4
	
	input:
	path pop_map
	tuple val(species), val(library_prep), path(vcf)
	
	output:
	path "*"
	
	when:
	params.variant_call_only == false && params.selection_scan == true
	
	script:
	"""
	ohana -geno input_file.geno \
	-pos input_file.pos \
	-pop ${pop_map} \
	-out output_file \
	-sel BAYENV \
	-threads ${task.cpus} \
	-minmaf ${params.minor_allele_frequency} \
	-mincoverage ${params.min_depth} \
	-n 5000
	"""
}


process ANGSD_GWAS {
	
	/* 
	This process does something described here
	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(bam_files)
	
	output:
	path "*.assoc"
	
	when:
	params.variant_call_only == false && params.genome_wide_association == true
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	angsd -bam input_list.txt \
	-out ${species}_${library_prep} \
	-doAsso 2 \
	-GL 1 \
	-doMaf 1 \
	-SNP_pval 1e-6 \
	-minMapQ 30 \
	-minQ 20 \
	-minInd ${params.max_snp_missingness} \
	-minMaf ${params.minor_allele_frequency} \
	-setMaxDepth ${params.max_depth}
	"""
}


process NGSLD {
	
	/* 
	This process does something described here
	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(vcf)
	
	output:
	path "*.ld"
	
	when:
	params.variant_call_only == false && params.linkage_disequilibrium == true
	
	script:
	"""
	NGSLD -vcf ${vcf} \
	-out "${species}_${library_prep}" \
	-ld \
	-maf ${params.minor_allele_frequency} \
	-ref ${params.reference}
	"""
}


process ROH {
	
	/* 
	This process does something described here
	*/
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(bam)
	
	output:
	path ".txt"
	
	when:
	params.esimate_roh == true && params.variant_call_only == false
	
	script:
	"""
	samtools view -H ${bam} | grep '^@RG' | cut -f 2 | cut -d ":" -f 2 > sample_list.txt && \
	samtools roh \
	-r chrom:start-end \
	-c sample_list.txt \
	-o ${sample}_roh.txt \
	input_file.bam
	"""
}


process PEDIGREE_STRUCTURES {
	
	/* 

	This process was inspired by the Molecular Ecology Resources study 
	Petty et al. 2019, titled:

	"Pedigree reconstruction and distant pairwise relatedness estimation 
	from genome sequence data: A demonstration in a population of rhesus 
	macaques (Macaca mulatta)"
	
	Code for this process comes from that study's GitHub repository, which 
	may be viewed at the link below:

	https://github.com/belowlab/tnprc-pedigrees
	
	*/
	
	tag "${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(vcf)
	
	output:
	path "*"
	
	when:
	params.pedigree_estimation == true && params.variant_call_only == false
	
	shell:
	"""
	pedigree_structures.sh ${vcf}
	"""
}


process D_STATISTIC {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	params.variant_call_only == false
	
	script:
	"""
	
	"""
}


// Population-level analyses

process ANGSD_AF {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:
	
	
	when:
	params.variant_call_only == false
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	angsd -b bam_list.txt \
	-doMajorMinor 1 \
	-doMaf 1 \
	-out ${species}_${library_prep} \
	-P ${task.cpus}
	"""
}


process ANGSD_SFS {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:
	path "*"
	
	when:
	params.variant_call_only == false && params.site_frequency_spectra == true
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	angsd -b bam_list.txt \
	-doSaf 1 \
	-ref ${params.reference} \
	-GL 1 \
	-nThreads ${task.cpus} \
    -out ${species}_${library_prep} && \
	realSFS ${species}_${library_prep}.saf.idx \
	-maxIter 100 \
	-tole 1e-6 -P ${task.cpus} > ${species}_${library_prep}_sfs.txt
	"""
}


process VISUALIZE_SFS {
	
	/* 
	This process does something described here
	*/
	
	tag "${prep} ${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple path(sfs), val(species), val(prep), val(sample_size)
	
	output:
	path "*.pdf"
	
	script:
	"""
	SFS_plotting.R ${sfs} ${species} ${params.date}
	"""
}


process STAIRWAY_PLOT {
	
	/* 
	This process does something described here
	*/
	
	tag "${prep} ${species}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple path(sfs), val(species), val(prep), val(sample_size)
	
	output:
	path "*${species}_${prep}*"
	
	when:
	params.variant_call_only == false && params.site_frequency_spectra == true && stairwayplot == true
	
	script:
	"""
	
	# create stairway plot blueprint based on VCF and settings in nextflow.config
	create_stairwayplot_blueprint.R ${sfs} \
	${species} \
	${species} \
	${prep} \
	${sample_size} \
	${params.genome_length} \
	${params.year_per_generation} \
	${params.mutation_rate} \
	${params.random_seed} \
	${params.whether_folded}
	
	# Pull stairway plot files
	wget https://github.com/xiaoming-liu/stairway-plot-v2/raw/master/stairway_plot_v2.1.1.zip && \
	unzip -o stairway_plot_v2.1.1.zip -d . && \
	rm stairway_plot_v2.1.1.zip && \
	mv stairway_plot_v2.1.1/stairway_plot_es/ . && \
	rm -rf stairway_plot_v2.1.1/
	
	# create stairway plot shell script and then run it
	java -cp stairway_plot_es Stairbuilder ${species}_${species}_${prep}.blueprint && \
	bash ${species}_${species}_${prep}.blueprint.sh
	
	# error out if summary files don't exist
	if [ \$(find . -maxdepth 1 -type f -name "*.final.summary*" | wc -l) -eq 0 ]; then
		error 1
	fi
	
	"""
}


process ANGSD_FST {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	params.variant_call_only == false
	
	script:
	"""
	
	"""
}


process NGSTOOLS_FST {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	params.variant_call_only == false
	
	script:
	"""
	
	"""
}


process VCFLIB_ASSESS_FST {
	
	/* 
	This process does something described here
	*/
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	params.variant_call_only == false
	
	script:
	"""
	
	"""
}


// PROCESS CODE TEMPLATE:
// process PROCESS_NAME {
	
// 	/* 
// 	This process does something described here
// 	*/
	
// 	tag "${tag}"
// 	publishDir params.results, mode: 'copy'
	
// 	memory 1.GB
// 	cpus 1
// 	time '10minutes'
	
// 	input:
	
	
// 	output:
	
	
// 	when:
	
	
// 	script:
// 	"""
	
// 	"""
// }

// --------------------------------------------------------------- //