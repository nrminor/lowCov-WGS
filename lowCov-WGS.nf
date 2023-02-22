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
    // Basic bioinformatic processing and QC
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
        FASTQC.out
    )

    MAP_TO_REFERENCE (
        TRIM_ADAPTERS.out
    )

    ASSESS_DEPTH (
        MAP_TO_REFERENCE.out
    )

    ANGSD_GL (
        MAP_TO_REFERENCE.out
    )

    // Individual-level analyses
    STRUCTURE ()

    NGSADMIX ()

    PCANGSD ()

    FASTPCA ()

    OHANA ()

    ANGSD_GWAS ()

    NGSLD ()

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

process MERGE_READS {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path(reads2_path)
	
	output:
    tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	
	script:
    if ( reads2_path.exists() ){
        """
        bbmerge.sh ${reads1_path} ${reads2_path}
        """
    } else {
        """
        mv ${reads1_path} "se_${reads1_path}"
        """
    }
	
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