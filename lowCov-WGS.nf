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
	
	
	// WORKFLOW STEPS
    // Routine bioinformatic processing and QC
    MERGE_READS (
        ch_reads
    )

    ORIENT_READS (
        MERGE_READS.out
    )

    FASTP_FILTER (
        ORIENT_READS.out
    )

    TRIM_ADAPTERS (
        FASTP_FILTER.out
    )

    FASTQC (
        TRIM_ADAPTERS.out
    )

    MULTIQC (
        FASTQC.out.collect()
    )

    MAP_TO_REFERENCE (
        TRIM_ADAPTERS.out
    )

    ASSESS_DEPTH (
        MAP_TO_REFERENCE.out
    )

    PLOT_DEPTH (
        ASSESS_DEPTH.out
    )

    ANGSD_GL (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> sample, population, species, library_prep, bam }
            .groupTuple( by:[2,3] )
    )

    CALL_VARIANTS (
        MAP_TO_REFERENCE.out
    )

    FILTER_VARIANTS (
        CALL_VARIANTS.out
    )

    MERGE_VARIANTS (
        FILTER_VARIANTS.out.vcf
			.map { species, library_prep, vcf -> tuple( file(raw_vcf), species, prep_type ) }
			.filter { file(it[0]).countLines() > 0 }
			.groupTuple( by: [1,2] ),
		FILTER_VARIANTS.out.index.collect()
    )

    // Individual-level analyses
    STRUCTURE (
        MERGE_VARIANTS.out
    )

    NGSADMIX (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
    )

    PCANGSD ()

    FASTPCA ()

    OHANA ()

    ANGSD_GWAS ()

    NGSLD ()

    ROH ()

    // D_STATISTIC ()

    // Population-level analyses
    ANGSD_AF ()
    
    ANGSD_SFS ()

    VISUALIZE_SFS ()

    BUILD_STAIRWAY_PLOT_SCRIPT ()

    STAIRWAY_PLOT ()

    ANGSD_FST ()

    NGSTOOLS_FST ()

    VCFLIB_ASSESS_FST ()
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //


// Basic bioinformatic processing and QC

process MERGE_READS {
	
	// This process does something described here
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path(reads2_path)
	
	output:
    tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	
	script:
    if ( reads2_path.exists() ){
        """
        bbmerge.sh in1=${reads1_path} in2=${reads2_path} \
        out=${sample}_${species}_${library_prep}.fastq.gz

        """
    } else {
        """
        mv ${reads1_path} ${sample}_${species}_${library_prep}.fastq.gz
        """
    }
	
}


process ORIENT_READS {
	
	// This process does something described here
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	input:
    tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	script:
	"""
    vsearch --orient ${reads} \
    --output ${sample}_${species}_${library_prep}_oriented.fastq.gz \
    --fastq_ascii 33
	"""
}


process FASTP_FILTER {
	
	// This process does something described here
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	script:
	"""
	fastp --in1 ${reads} \
    --out1 ${sample}_${species}_${library_prep}_filtered.fastq.gz \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --detect_adapter_for_pe --detect_adapter --correction \
    --trim_tail1 5 --trim_tail2 5 \
    --thread ${task.cpus}
	"""
}


process TRIM_ADAPTERS {
	
	// This process does something described here
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	script:
	"""
	java -jar trimmomatic-0.35.jar SE -phred33 ${reads} \
    ${sample}_${species}_${library_prep}_trimmed.fastq.gz \
    ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50
	"""
}


process FASTQC {
	
	// This process does something described here
	
	tag "${sample}"
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	path "*"
	
	script:
	"""
	fastqc -f fastq -o ./ -t ${task.cpus} --extract ${reads}
	"""
}


process MULTIQC {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
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
	
	// This process does something described here
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	cpus 4
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.bam")
	
	script:
	"""
	bwa mem -t ${task.cpus} ${params.reference} ${reads} \
    | samtools view -bS - > ${sample}_${species}_${library_prep}.bam
	"""
}


process ASSESS_DEPTH {
	
	// This process does something described here
	
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(bam)
	
	output:
	path "*"
	
	script:
	"""
	mosdepth -t ${task.cpus} ${sample}_${species}_${library_prep} ${bam}
	"""
}


process PLOT_DEPTH {
	
	// This process does something described here
	
	publishDir params.results, mode: 'copy'
	
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
	
	// This process does something described here
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
	cpus 8
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:
	
	
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
	
	// This process does something described here
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(samples), val(populations), val(species), val(library_prep), path(bam_files)
	
	output:
	tuple val(species), val(library_prep), path("*.vcf.gz")
	
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


process FILTER_VARIANTS {
	
	// This process does something described here
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(vcf)
	
	output:
	tuple val(species), val(library_prep), path("*.vcf.gz"), emit: vcf
	path "*.tbi", emit: index
	
	script:
	"""
    
	"""
}


process MERGE_VARIANTS {
	
	// This process does something described here
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
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
	
	sample_ids=`bcftools query -l !{species}_multisample_!{prep}.vcf.gz`
	sample_size=`bcftools query -l !{species}_multisample_!{prep}.vcf.gz | wc -l`
	minInd="$((${sample_size} - !{params.max_snp_missingness}))"
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / ${sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: ${sample_ids}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	
	'''
}




// Individual-level analyses

process STRUCTURE {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple path(merged_vcf), val(species), val(library_prep), val(sample_size)
	
	output:
	
	
	when:
	params.structure == true
	
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
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(species), val(library_prep), path(bam_files)
	
	output:
	
	
	script:
	"""
    find . -name "*.bam" > bam_file_list.txt && \
	ngsAdmix -P bam_file_list.txt \
    -B ${params.reference} \
    -K 3 \
    -o ${species}_${library_prep}
	"""
}


process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}


process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

// --------------------------------------------------------------- //