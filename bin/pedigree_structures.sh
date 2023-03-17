#!/bin/bash

input_vcf=$1

gatk VariantFiltration \
--variant ${input_vcf} \
--output allchr_filtered.vcf \
--filter-name "hard_filter" \
--filter-expression "AC == 0 || QD < 10.0 || FS > 2.0 || ( MQ < 59.0 && MQ > 61.00 ) || SOR > 1.5"

# converting vcf to PLINK format files for PRIMUS:
plink --vcf allchr_filtered.vcf \
    --make-bed --out allchr_filtered_plink

# running PRIMUS:
run_PRIMUS.pl --file allchr_filtered_plink \
    --genome --keep_inter_files --internal_ref \
    --degree_rel_cutoff 1 --output_dir PRIMUS_results

# combining networks in PRIMUS:
for network in {0..29}
do
    grep -P "^${network}\t" > group${network}_keep
    plink --bfile allchr_filtered_plink \
        --keep group${network}_keep --make-bed \
        --out allchr_filtered_plink_group${network}
    run_PRIMUS.pl --file allchr_filtered_plink_group${network} \
        --genome --keep_inter_files --degree_rel_cutoff 2 \
        --internal_ref --output_dir PRIMUS_group${network}
done

# splitting by chromosome:
for c in {1..20}
do
    plink --bfile allchr_filtered_plink \
        --chr ${c} --recode vcf \
        --out filtered_chr${c}
done

# phasing:
for c in {1..20}
do
    echo "java -jar beagle.25Nov19.28d.jar gt=filtered_chr${c}.vcf nthreads=3 out=filtered_chr${c}_phased"
done | parallel

# converting phased to PLINK format, preparing GERMLINE input:
for c in {1..20}
do
    vcftools --gzvcf filtered_chr${c}_phased.vcf.gz \
        --plink --out filtered_chr${c}_phased_plink
    awk '{ OFS="\t"; $3=($4/1000000)*0.433; print; }' filtered_chr${c}_phased_plink.map > filtered_chr${c}_phased_plink_cm.map
    echo "1" > germline_chr${c}.run
    echo "filtered_chr${c}_phased_plink_cm.map" >> germline_chr${c}.run
    echo "filtered_chr${c}_phased_plink.ped" >> germline_chr${c}.run
    echo "filtered_chr${c}_phased_plink_germline" >> germline_chr${c}.run
done

# running GERMLINE:
het=1
hom=2
for c in {1..20}
do
    echo "germline -min_m 2.5 -err_het ${het} -err_hom ${hom} < germline_chr${c}.run"
done | parallel

# running ERSA
ersa --segment_files=*.match --number_of_chromosomes=20 \
    --rec_per_meioses=13.6239623 --confidence_level=0.999 \
    --output_file=ersa_results.txt --model_output_file=ersa_models.txt

# varying recombinations per meiosis in ERSA:
for r in $(seq 6.0 0.2 13.4)
do
        ersa --segment_files=*.match --number_of_chromosomes=20 \
        --confidence_level=0.999 --control_files=non0419.match \
        --mask_common_shared_regions=true --mask_region_threshold=6 \
        --rec_per_meioses=${r} --output_file=ersa_results_varied_rpm_${r}.txt \
        --model_output_file=ersa_models_varied_rpm_${r}.txt
done

# running PADRE
run_PRIMUS.pl --project_summary PRIMUS_results/Summary_allchr_filtered_plink.genome.txt \
    --ersa_model_output ersa_models.txt --ersa_results ersa_results.txt --degree_rel_cutoff 1


## Downsampled analyses
# variant filtering for downsampled data
bcftools filter -i 'FORMAT/DP>5 && GQ>12' -S . \
    -o downsampled_5.0_filt.vcf starting_vcf_file.vcf

# converting vcf to PLINK format files for PRIMUS:
plink --vcf downsampled_5.0_filt.vcf \
    --make-bed --out downsampled_5.0_filt_plink

# running PRIMUS:
run_PRIMUS.pl --file downsampled_5.0_filt_plink \
    --genome --keep_inter_files --internal_ref \
    --degree_rel_cutoff 1 --output_dir PRIMUS_results_5.0

# combining networks in PRIMUS:
for network in {0..29}
do
    plink --bfile downsampled_5.0_filt_plink \
        --keep group${network}_keep --make-bed \
        --out downsampled_5.0_filt_plink_group${network}
    run_PRIMUS.pl --file downsampled_5.0_filt_plink_group${network} \
        --genome --keep_inter_files --degree_rel_cutoff 2 \
        --internal_ref --output_dir PRIMUS_results_5.0_group${network}
done

# splitting by chromosome:
for c in {1..20}
do
    plink --bfile downsampled_5.0_filt_plink \
        --chr ${c} --recode vcf \
        --out downsampled_5.0_filt_chr${c}
done

# phasing:
for c in {1..20}
do
    echo "java -jar beagle.25Nov19.28d.jar gt=downsampled_5.0_filt_chr${c}.vcf nthreads=3 out=downsampled_5.0_filt_chr${c}_phased"
done | parallel

# converting phased to PLINK format, preparing GERMLINE input:
for c in {1..20}
do
    vcftools --gzvcf downsampled_5.0_filt_chr${c}_phased.vcf.gz \
        --plink --out downsampled_5.0_filt_chr${c}_phased_plink
    awk '{ OFS="\t"; $3=($4/1000000)*0.433; print; }' downsampled_5.0_filt_chr${c}_phased_plink.map > downsampled_5.0_filt_chr${c}_phased_plink_cm.map
    echo "1" > germline_chr${c}.run
    echo "downsampled_5.0_filt_chr${c}_phased_plink_cm.map" >> germline_chr${c}.run
    echo "downsampled_5.0_filt_chr${c}_phased_plink.ped" >> germline_chr${c}.run
    echo "downsampled_5.0_filt_chr${c}_phased_germline" >> germline_chr${c}.run
done

# running GERMLINE:
het=19
hom=5
for c in {1..20}
do
    echo "germline -min_m 2.5 -err_het ${het} -err_hom ${hom} < germline_chr${c}.run"
done | parallel

# running ERSA
ersa --segment_files=*.match --number_of_chromosomes=20 \
    --rec_per_meioses=13.6239623 --confidence_level=0.999 \
    --output_file=ersa_results_5.0.txt --model_output_file=ersa_models_5.0.txt

# running PADRE
run_PRIMUS.pl --project_summary PRIMUS_results_5.0/Summary_downsampled_5.0_filt_plink.genome.txt \
    --ersa_model_output ersa_models_5.0.txt --ersa_results ersa_results_5.0.txt --degree_rel_cutoff 1