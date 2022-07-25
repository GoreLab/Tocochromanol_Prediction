## filter DR2>0.8
cp /workdir/dw524/Hapmap3_Ames_Gorelab/vitamaize/Ames_vitamaize_imputation_stats.txt .
awk '{ if($4 >= 0.8) { print }}' Ames_vitamaize_imputation_stats.txt > Ames_imputation_stats_DR2.txt

n_chromosomes=10
prefix=AGPv4_Ames_vitamaize

## filter for bi-allelic, DR2, and only chosen taxa
mkdir 1.DR2
echo 'step 1...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf /workdir/dw524/Hapmap3_Ames_Gorelab/vitamaize/${prefix}_chr${i}.imputed.vcf.gz \
--positions Ames_imputation_stats_DR2.txt \
--min-alleles 2 --max-alleles 2 --recode --stdout \
| gzip -c > ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz &
done

wait

###filter maf, LD pruning
mkdir 3.filtering
echo 'filter maf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz --maf 0.01 --recode --stdout | gzip -c > ./3.filtering/${prefix}${i}_maf.vcf.gz &
done 

wait

# echo 'filter LD...'

# for i in $(seq 1 1 ${n_chromosomes}) 
# do 
# /programs/plink-1.9-x86_64-beta5/plink --vcf ./3.filtering/${prefix}${i}_maf.vcf.gz --indep-pairwise 100 25 0.9 --out ./3.filtering/${prefix}${i}_maf_LD &
# done 

# wait

# echo 'rewrite vcf...'
# for i in $(seq 1 1 ${n_chromosomes}) 
# do 
# vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --snps ./3.filtering/${prefix}${i}_mac_LD.prune.in --recode --stdout | gzip -c > ./3.filtering/${prefix}${i}_maf_LD.vcf.gz &
# done 

# wait

## convert to hmp.txt for GAPIT 
mkdir 4.hmp
echo 'vcf to hmp...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -fork1 \
-vcf ./3.filtering/${prefix}${i}_maf.vcf.gz -export ./4.hmp/${prefix}${i}_maf_full \
-exportType Hapmap 
done 

wait 

###get summary statistics before
mkdir 2.summary
echo 'summary...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g \
-vcf ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz -genotypeSummary taxa,site \
-export ./2.summary/${prefix}${i}_full
done 

wait

###get summary statistics after
mkdir 5.summary

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g \
-vcf ./3.filtering/${prefix}${i}_maf.vcf.gz -genotypeSummary taxa,site \
-export ./5.summary/${prefix}${i}_full_after
done 


##convert to numerical
mkdir 6.numerical
echo 'vcf to numerical...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --012 --out 6.numerical/${prefix}${i}_maf &
done



##LD pruning for MLMM
mkdir 7.LD_pruning_for_MLMM
echo 'filter LD...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf 3.filtering/${prefix}${i}_maf.vcf.gz --indep-pairwise 100 25 0.99 --out ./7.LD_pruning_for_MLMM/${prefix}${i}_maf_LD &
done 

wait

echo 'rewrite vcf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --snps ./7.LD_pruning_for_MLMM/${prefix}${i}_maf_LD.prune.in --recode --stdout | gzip -c > ./7.LD_pruning_for_MLMM/${prefix}${i}_maf_LD.vcf.gz &
done 

wait

echo 'convert to 012'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./7.LD_pruning_for_MLMM/${prefix}${i}_maf_LD.vcf.gz --012 --out ./7.LD_pruning_for_MLMM/${prefix}${i}_maf_LD &
done 


##get SNPs for kinship
mkdir 8.kinship

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf 3.filtering/${prefix}${i}_maf.vcf.gz --indep-pairwise 100 25 0.2 --out ./8.kinship/${prefix}${i}_maf_LD0.2 &
/programs/plink-1.9-x86_64-beta5/plink --vcf 3.filtering/${prefix}${i}_maf.vcf.gz --indep-pairwise 100 25 0.1 --out ./8.kinship//${prefix}${i}_maf_LD0.1 &
done 
wait


for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --snps ./8.kinship/${prefix}${i}_maf_LD0.2.prune.in --recode --out ./8.kinship/${prefix}${i}_maf_LD0.2 &
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --snps ./8.kinship/${prefix}${i}_maf_LD0.1.prune.in --recode --out ./8.kinship/${prefix}${i}_maf_LD0.1 &
done
wait

for i in $(seq 1 1 ${n_chromosomes})
do
bgzip -f ./8.kinship/${prefix}${i}_maf_LD0.1.recode.vcf &
bgzip -f ./8.kinship/${prefix}${i}_maf_LD0.2.recode.vcf &
done
wait

for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf ./8.kinship/${prefix}${i}_maf_LD0.1.recode.vcf.gz &
tabix -f -p vcf ./8.kinship/${prefix}${i}_maf_LD0.2.recode.vcf.gz &
done
wait

#merge
bcftools concat ./8.kinship/${prefix}*_maf_LD0.1.recode.vcf.gz -Oz -o ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz &
bcftools concat ./8.kinship/${prefix}*_maf_LD0.2.recode.vcf.gz -Oz -o ./8.kinship/${prefix}_LD0.2_all_chr.vcf.gz &
wait
tabix -f -p vcf ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz &
tabix -f -p vcf ./8.kinship/${prefix}_LD0.2_all_chr.vcf.gz &

vcftools --gzvcf ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz --012 --out ./8.kinship/${prefix}_LD0.1_all_chr &
vcftools --gzvcf ./8.kinship/${prefix}_LD0.2_all_chr.vcf.gz --012 --out ./8.kinship/${prefix}_LD0.2_all_chr &



##LD pruning for eQTL
mkdir 9.for_eQTL
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz \
--keep expression_sample_list_for_eqtl_first_run.txt --min-alleles 2 --max-alleles 2 \
--maf 0.01 --recode --out 9.for_eQTL/${prefix}${i}_exp_subset &
done

wait 
 
for i in $(seq 1 1 ${n_chromosomes}) 
do 
bgzip -f ./9.for_eQTL/${prefix}${i}_exp_subset.recode.vcf &
done

echo 'filter LD...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf ./9.for_eQTL/${prefix}${i}_exp_subset.recode.vcf.gz --indep-pairwise 100 25 0.8 --out ./9.for_eQTL/${prefix}${i}_maf_LD0.8 &
/programs/plink-1.9-x86_64-beta5/plink --vcf ./9.for_eQTL/${prefix}${i}_exp_subset.recode.vcf.gz --indep-pairwise 100 25 0.7 --out ./9.for_eQTL/${prefix}${i}_maf_LD0.7 &
done 

wait

echo 'rewrite vcf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./9.for_eQTL/${prefix}${i}_exp_subset.recode.vcf.gz --snps ./9.for_eQTL/${prefix}${i}_maf_LD0.8.prune.in --recode --out ./9.for_eQTL/${prefix}${i}_maf_LD0.8 &
vcftools --gzvcf ./9.for_eQTL/${prefix}${i}_exp_subset.recode.vcf.gz --snps ./9.for_eQTL/${prefix}${i}_maf_LD0.7.prune.in --recode --out ./9.for_eQTL/${prefix}${i}_maf_LD0.7 &
done 

wait

for i in $(seq 1 1 ${n_chromosomes}) 
do 
bgzip -f ./9.for_eQTL/${prefix}${i}_maf_LD0.8.recode.vcf &
bgzip -f ./9.for_eQTL/${prefix}${i}_maf_LD0.7.recode.vcf &
done

wait

for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf ./9.for_eQTL/${prefix}${i}_maf_LD0.8.recode.vcf.gz &
tabix -f -p vcf ./9.for_eQTL/${prefix}${i}_maf_LD0.7.recode.vcf.gz &
done
wait

#merge
bcftools concat ./9.for_eQTL/${prefix}*_maf_LD0.8.recode.vcf.gz -Oz -o ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8.recode.vcf.gz &
bcftools concat ./9.for_eQTL/${prefix}*_maf_LD0.7.recode.vcf.gz -Oz -o ./9.for_eQTL/${prefix}_all_chr_maf_LD0.7.recode.vcf.gz 

#index
tabix -f -p vcf ./9.for_eQTL/${prefix}_all_chr_maf_LD0.7.recode.vcf.gz &
tabix -f -p vcf ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8.recode.vcf.gz 

echo 'convert to hmp'
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -fork1 \
-vcf ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8.recode.vcf.gz -export ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8 \
-exportType Hapmap &

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -fork1 \
-vcf ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8.recode.vcf.gz -export ./9.for_eQTL/${prefix}_all_chr_maf_LD0.8 \
-exportType Hapmap &



mkdir 10.kinship_for_TWAS_v1_1
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz \
--keep line_names_for_TWAS_v1.1.txt --min-alleles 2 --max-alleles 2 \
--maf 0.01 --recode --out 10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset &
done

wait 
 
for i in $(seq 1 1 ${n_chromosomes}) 
do 
bgzip -f 10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf &
done

echo 'change to 012'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf 10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf.gz \
 --012 --out 10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset &
done

wait

echo 'filter LD...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf.gz --indep-pairwise 100 25 0.1 --out ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.1 &
done 

wait

echo 'rewrite vcf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf.gz --snps ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.1.prune.in --recode --out ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.1 &
done 

wait

for i in $(seq 1 1 ${n_chromosomes}) 
do 
bgzip -f ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.1.recode.vcf &
done

wait

for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.1.recode.vcf.gz &
done
wait

#merge
bcftools concat ./10.kinship_for_TWAS_v1_1/${prefix}*_exp_maf_LD0.1.recode.vcf.gz -a -Oz -o ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1.recode.vcf.gz

wait

#index
tabix -f -p vcf ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1.recode.vcf.gz &

wait 

echo 'convert to hmp'
/programs/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -Xmx300g \
-inputFile ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1.recode.vcf.gz \
-outputFile ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1.sorted

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -fork1 \
-vcf ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1.sorted.vcf -export ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.1 \
-exportType Hapmap &

#directly subset from kinship set
vcftools --gzvcf 8.kinship/AGPv4_Ames_vitamaize_LD0.1_all_chr.vcf.gz --keep line_names_for_TWAS_v1.1.txt --maf 0.01 --min-alleles 2 --max-alleles 2 --recode --out AGPv4_Ames_vitamaize_LD0.1_all_chr_TWAS_eQTL_subset


# 0.99 filtering
for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf.gz --indep-pairwise 100 25 0.99 --out ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.99 &
done 

wait

echo 'rewrite vcf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_subset.recode.vcf.gz --snps ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.99.prune.in --recode --out ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.99 &
done 

wait

for i in $(seq 1 1 ${n_chromosomes}) 
do 
bgzip -f ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.99.recode.vcf &
done

wait

for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf ./10.kinship_for_TWAS_v1_1/${prefix}${i}_exp_maf_LD0.99.recode.vcf.gz &
done
wait

#merge
bcftools concat ./10.kinship_for_TWAS_v1_1/${prefix}*_exp_maf_LD0.99.recode.vcf.gz -a -Oz -o ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.99.recode.vcf.gz

wait

#index
tabix -f -p vcf ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.99.recode.vcf.gz &

wait 

echo 'convert to numerical'
vcftools --gzvcf ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.99.recode.vcf.gz \
--012 --out ./10.kinship_for_TWAS_v1_1/${prefix}_all_chr_exp_maf_LD0.99 &




