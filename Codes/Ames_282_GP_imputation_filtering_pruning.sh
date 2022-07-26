##download raw SNP info
iget -K /iplant/home/shared/panzea/genotypes/GBS/v27/ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz

###################################################################
# 1. get corresponding raw genotype and filter for biallelic
###################################################################
##Ames
cp ../../../vitamaize/genotype/ames_raw_1462_443419_consensus.recode.vcf . &

##282
vcftools --gzvcf ../../../ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz --keep 282.taxa.list.final.2020.02.06.txt \
--max-alleles 2 --min-alleles 2 --maf 0.000001 \
--recode --out 282_unimputed_unfiltered &

wait

###################################################################
# 2. get overlapping SNPs and file format
###################################################################
##
cut -f 3 ames_raw_1462_443419_consensus.recode.vcf > Ames_pos.txt
cut -f 3 282_unimputed_unfiltered.recode.vcf > 282_pos.txt

grep -Fxf "282_pos.txt" "Ames_pos.txt" > comm_pos.txt
wc -l comm_pos.txt

#bgzip and index
bgzip -f ames_raw_1462_443419_consensus.recode.vcf &
bgzip -f 282_unimputed_unfiltered.recode.vcf &
wait

tabix -f -p vcf ames_raw_1462_443419_consensus.recode.vcf.gz &
tabix -f -p vcf 282_unimputed_unfiltered.recode.vcf.gz

vcftools --gzvcf ames_raw_1462_443419_consensus.recode.vcf.gz --snps comm_pos.txt --recode --out Ames_unimputed_unfiltered_overlap &
vcftools --gzvcf 282_unimputed_unfiltered.recode.vcf.gz --snps comm_pos.txt --recode --out 282_unimputed_unfiltered_overlap &

wait


#change to hmp
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf Ames_unimputed_unfiltered_overlap.recode.vcf \
-export -exportType Hapmap &

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -vcf 282_unimputed_unfiltered_overlap.recode.vcf \
-export -exportType Hapmap &

wait

###################################################################
# 3. merge genotype files
###################################################################
#using R script to merge two files
Rscript merging_Hapmap_files.R

#change to vcf
/programs/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -Xmx300g \
-inputFile merged_ames_282_biallelic.hmp.txt -outputFile merged_ames_282_biallelic_sorted.hmp.txt 
rm merged_ames_282_biallelic.hmp.txt

/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -h merged_ames_282_biallelic_sorted.hmp.txt \
-export -exportType VCF

#removing het
sed -i '/^#/! s/0\/1/\.\/\./g; s/1\/0/\.\/\./g' merged_ames_282_biallelic_sorted.vcf

#bgzip and index
bgzip -f merged_ames_282_biallelic_sorted.vcf 
tabix -f -p vcf merged_ames_282_biallelic_sorted.vcf.gz

#normalizing
bcftools norm -m +any -Oz -o merged_ames_282_biallelic_sorted-normalized.vcf.gz merged_ames_282_biallelic_sorted.vcf.gz 

#filtering
bcftools view -m2 -M2 -v snps -Oz -o merged_ames_282_biallelic-unsorted.vcf.gz merged_ames_282_biallelic_sorted-normalized.vcf.gz
rm merged_ames_282_biallelic_sorted-normalized.vcf.gz

# Sorting
bcftools sort -Oz -o merged_ames_282_biallelic.vcf.gz merged_ames_282_biallelic-unsorted.vcf.gz 
rm merged_ames_282_biallelic-unsorted.vcf.gz

###################################################################
# 4. Separating chromosomes and indexing VCF
###################################################################
cd imputation

beagle_version=/workdir/dw524/programs/beagle.28Sep18.793.jar
n_chromosomes=10

echo Indexing...

#vitamaize
for i in $(seq 1 1 ${n_chromosomes})
do
# Converting
vcftools --gzvcf ../merging/merged_ames_282_biallelic.vcf.gz \
--chr ${i} --recode --out merged_ames_282_chr${i} &
done

wait 

# zipping
for i in $(seq 1 1 ${n_chromosomes})
do
bgzip -f merged_ames_282_chr${i}.recode.vcf &
done
wait

#indexing
for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf merged_ames_282_chr${i}.recode.vcf.gz &
done



###################################################################
# 5. imputation
###################################################################
input_dir=/workdir/dw524/Hapmap3_Ames_Gorelab/
output_dir=/workdir/dw524/tocochromanol_GP/imputation/v2/imputation/
Ames_prefix=${output_dir}merged_ames_282_chr
ref_prefix=${input_dir}Hmp321/hmp321_282_agpv4_merged_chr

echo Imputing in Ames...
for i in $(seq 1 1 ${n_chromosomes})
do
echo Chromosome ${i}
ref=${ref_prefix}${i}.imputed
study=${Ames_prefix}${i}.recode
out=${Ames_prefix}${i}.imputed
java -Xmx200g -jar ${beagle_version} ref=${ref}.vcf.gz gt=${study}.vcf.gz burnin=10 \
nthreads=80 map=${input_dir}NAM_genetic_map/Map_NAM_Chr${i}_cleaned.txt iterations=15 \
ne=50000 out=${out} 
done

#indexing
for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf merged_ames_282_chr${i}.imputed.vcf.gz &
done

wait

#imputation accuracy
echo -e "CHROM\tPOS\tAF_Ames\tR2_Beagle" > merged_ames_282_imputation_stats.txt

for i in $(seq 1 1 ${n_chromosomes})
do
bcftools query -f '%CHROM\t%POS\t%AF\t%DR2\n' merged_ames_282_chr${i}.imputed.vcf.gz >> merged_ames_282_imputation_stats.txt
done

###################################################################
# 6. filtering
###################################################################
cd ../filtering
## filter DR2>0.8
awk '{ if($4 >= 0.8) { print }}' ../imputation/merged_ames_282_imputation_stats.txt > merged_ames_282_imputation_stats_DR2.txt

n_chromosomes=10
prefix=merged_ames_282

## filter for bi-allelic, DR2, and only chosen taxa
mkdir 1.DR2
echo 'step 1, DR2...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ../imputation/${prefix}_chr${i}.imputed.vcf.gz \
--positions merged_ames_282_imputation_stats_DR2.txt \
--min-alleles 2 --max-alleles 2 --recode --stdout \
| gzip -c > ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz &
done

wait

###filter maf
mkdir 3.filtering
echo 'filter maf...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz --maf 0.01 --recode --stdout | gzip -c > ./3.filtering/${prefix}${i}_maf.vcf.gz &
done 

wait

## convert to hmp.txt 
mkdir 4.hmp
echo 'vcf to hmp...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g -fork1 \
-vcf ./3.filtering/${prefix}${i}_maf.vcf.gz -export ./4.hmp/${prefix}${i}_maf_full \
-exportType Hapmap &
done 

wait 

###get summary statistics before
mkdir 2.summary
echo 'summary...'
for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g \
-vcf ./1.DR2/${prefix}_chr${i}.imputed.full.vcf.gz -genotypeSummary taxa,site \
-export ./2.summary/${prefix}${i}_full &
done 

wait

###get summary statistics after
mkdir 5.summary

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/tassel-5-standalone/run_pipeline.pl -Xmx300g \
-vcf ./3.filtering/${prefix}${i}_maf.vcf.gz -genotypeSummary taxa,site \
-export ./5.summary/${prefix}${i}_full_after &
done 
wait

##convert to numerical
mkdir 6.numerical
echo 'vcf to numerical...'

for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --012 --out 6.numerical/${prefix}${i}_maf &
done


##LD pruning
mkdir 8.kinship

for i in $(seq 1 1 ${n_chromosomes}) 
do 
/programs/plink-1.9-x86_64-beta5/plink --vcf 3.filtering/${prefix}${i}_maf.vcf.gz --indep-pairwise 100 25 0.1 --out ./8.kinship//${prefix}${i}_maf_LD0.1 &
done 
wait


for i in $(seq 1 1 ${n_chromosomes}) 
do 
vcftools --gzvcf ./3.filtering/${prefix}${i}_maf.vcf.gz --snps ./8.kinship/${prefix}${i}_maf_LD0.1.prune.in --recode --out ./8.kinship/${prefix}${i}_maf_LD0.1 &
done
wait

for i in $(seq 1 1 ${n_chromosomes})
do
bgzip -f ./8.kinship/${prefix}${i}_maf_LD0.1.recode.vcf &
done
wait

for i in $(seq 1 1 ${n_chromosomes})
do
tabix -f -p vcf ./8.kinship/${prefix}${i}_maf_LD0.1.recode.vcf.gz &
done
wait

#merge
bcftools concat ./8.kinship/${prefix}*_maf_LD0.1.recode.vcf.gz -Oz -o ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz &
wait
tabix -f -p vcf ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz &

vcftools --gzvcf ./8.kinship/${prefix}_LD0.1_all_chr.vcf.gz --012 --out ./8.kinship/${prefix}_LD0.1_all_chr &





