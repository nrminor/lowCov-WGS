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
	REMOVE_OPTICAL_DUPLICATES (
		ch_reads
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

	QUALITY_TRIM (
		NORMALIZE_READS.out
	)

    MERGE_READS (
		QUALITY_TRIM.out,
		ch_reads
	)

    ORIENT_READS (
        MERGE_READS.out
    )

    FASTP_FILTER (
        ORIENT_READS.out
    )

    FASTQC (
        FASTP_FILTER.out
    )

    MULTIQC (
        FASTQC.out.collect()
    )

    MAP_TO_REFERENCE (
        FASTP_FILTER.out
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
        MERGE_VARIANTS.out.vcf
    )

    NGSADMIX (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
    )

    PCANGSD (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
	)

    FASTPCA (
		ANGSD_GL.out
	)

    OHANA (
        MERGE_VARIANTS.out.vcf
	)

    ANGSD_GWAS (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
	)

    NGSLD (
		FILTER_VARIANTS.out.vcf
	)

    ROH (
		MAP_TO_REFERENCE.out
	)

	PEDIGREE_STRUCTURES (
		MERGE_VARIANTS.out.vcf
	)

    // D_STATISTIC ()

    // Population-level analyses
    ANGSD_AF (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
	)
    
    ANGSD_SFS (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
	)

    VISUALIZE_SFS (
		ANGSD_SFS.out.sfs
	)

    STAIRWAY_PLOT (
		ANGSD_SFS.out.sfs
	)

    ANGSD_FST (
        MAP_TO_REFERENCE.out
            .map { sample, population, species, library_prep, seq_platform, bam -> species, library_prep, bam }
            .groupTuple( by:[0,1] )
	)

    NGSTOOLS_FST ()

    VCFLIB_ASSESS_FST (
		ANGSD_FST.out
			.mix (
				NGSTOOLS_FST.out
			)
	)
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

// Standard bioinformatic processing and QC. 
// -----------------------------------------
// These steps are designed for Illumina paired-end short reads and not to any other 
// platform or read length, but support for these sequence read configurations 
// may be added in the future.

process REMOVE_OPTICAL_DUPLICATES {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path(reads2_path)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	clumpify.sh in=${reads1_path} in2=${reads2_path} \
	out=${sample}_clumped.fastq.gz \
	threads=${task.cpus} \
	dedupe optical tossbrokenreads
	"""

}


process REMOVE_LOW_QUALITY_REGIONS {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	filterbytile.sh in=${reads} \
	out=${sample}_filtered_by_tile.fastq.gz \
	threads=${task.cpus}
	"""

}


process TRIM_ADAPTERS {
	
	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	bbduk.sh in=${reads} \
	out=${sample}_trim_adapters.fastq.gz \
	ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered \
	threads=${task.cpus}
	"""

}


process REMOVE_ARTIFACTS {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	bbduk.sh in=${reads} \
	out=${sample}_trim_adapters.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_ONE {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	bbmerge.sh in=${reads} \
	out=${sample}_remove_artifacts.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample}_ihist_merge1.txt \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_TWO {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	clumpify.sh in=${reads} \
	out=${sample}_eccc.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus}
	"""

}


process ERROR_CORRECT_PHASE_THREE {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	tadpole.sh in=${reads} \
	out=${sample}_ecct.fastq.gz \
	ecc k=62 ordered \
	threads=${task.cpus}
	"""

}


process NORMALIZE_READS {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	bbnorm.sh in=${reads} \
	out=${sample}_normalized.fastq.gz \
	target=100 \
	hist=${sample}_khist.txt \
	peaks=${sample}_peaks.txt
	"""

}


process MERGE_READS {
	
	/* 
	This process does something described here
	*/

	tag "${sample}"
	
	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	tuple val(no_QC_sample), val(no_QC_population), val(no_QC_species), val(no_QC_library_prep), val(no_QC_seq_platform), path(reads1_path), path(reads2_path)
	
	output:
    tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	
	script:
    if ( reads2_path.exists() ){
        """
        bbmerge-auto.sh in=${reads} \
		out=${sample}_${species}_${library_prep}.fastq.gz \
		outu=${sample}_unmerged.fastq.gz \
		strict k=93 extend2=80 rem ordered \
		ihist=${sample}_ihist_merge.txt \
		threads=${task.cpus}
        """
    } else {
        """
        mv ${reads} ${sample}_${species}_${library_prep}.fastq.gz
        """
    }
	
}


process QUALITY_TRIM {

	/* 
	This process does something described here
	*/

	tag "${sample}"

	cpus 8

	input:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)

	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads1_path), path("*.fastq.gz")

	script:
	"""
	bbduk.sh in=${reads} \
	out=${sample}_qtrimmed.fastq.gz \
	qtrim=r trimq=10 minlen=70 ordered \
	threads=${task.cpus}
	"""

}


process ORIENT_READS {
	
	/* 
	This process does something described here
	*/
	
	tag "${sample}"
	publishDir params.results, mode: 'copy'
	
	input:
    tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path(reads)
	
	output:
	tuple val(sample), val(population), val(species), val(library_prep), val(seq_platform), path("*.fastq.gz")
	
	script:
	"""
    vsearch --orient ${reads} \
	--db ${params.reference} \
    --output ${sample}_${species}_${library_prep}_oriented.fastq.gz \
    --fastq_ascii 33
	"""
}


process FASTP_FILTER {
	
	/* 
	This process does something described here
	*/
	
	tag "${sample}"
	
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


process FASTQC {
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
	tag "${species} ${library_prep}"
	publishDir params.results, mode: 'copy'
	
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
	
	/* 
	This process does something described here
	*/
	
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
	
	/* 
	This process does something described here
	*/
	
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