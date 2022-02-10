# Make Ames phenotype Data

# source
library(gdata)
library(data.table)

# define object
trait.all <- c("a.T", "d.T", "g.T",
               "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# file input
file.pheno.all <- "RAWDATA/Data_Ames/vitamaize_accession_removal_20191113.xlsx"
file.indv <- "RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv"
file.lambda <- "RAWDATA/Data_Ames/boxcox_transformation_applied_input_test_with_checks.txt"

# file output
dir.out <- "RESULT/1.1-MakeAmesPhenoData"
out.pheno <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv"
out.pheno.trans <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno_trans.csv"
out.info <- "RESULT/1.1-MakeAmesPhenoData/AmesInfo.csv"
out.BoxCoxParam <- "RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv"

# mkdir
dir.create(dir.out, recursive = TRUE)

# load phenotype info file
pheno.all <- read.xls(file.pheno.all, stringsAsFactors = F)


table(pheno.all$Removal_decision)

# phenotype data to be kept
tf <- (pheno.all$removed_vita == "no") & (pheno.all$Removal_decision == "keep")
table(tf) # TRUE should be 1462 -- OK!!

# 1462 phenotype data
df.info <- data.frame("ID" = pheno.all[tf, "Accession.Number_ion"],
                      "Pedigree" = pheno.all[tf, "Pedigree_vita"],
                      "Subpop.ISU.GB" = pheno.all[tf, "Comments_isu_gb_vita"],
                      "Seed.Weight" = pheno.all[tf, "seedweight"])
df.toco.pheno <- pheno.all[tf, c("Accession.Number_ion", trait.all)]
colnames(df.toco.pheno)[1] <- "ID"

# load ID of the genotype data
indv.tmp <- fread(file.indv, header = FALSE)
id.geno <- indv.tmp$V1

# check if all phenotyped accessions are in the genotye data
setdiff(df.toco.pheno$ID, id.geno) # perfect!

# load lambda
Lambda.all <- read.delim(file.lambda, stringsAsFactors = F)
lambda.vec <- Lambda.all[match(trait.all, Lambda.all$Trait), "Lambda"]
names(lambda.vec) <- trait.all

# make transformed phenotype data
const.vec <- lambda.vec # set a object
df.toco.pheno.trans <- df.toco.pheno # set a object
for (trait in trait.all) {
   # get values
   lambda <- lambda.vec[trait]
   pheno.vec <- df.toco.pheno[[trait]]
   
   # add constant so that al values are positive
   min.value <- min(pheno.vec, na.rm = T)
   if ( min.value < 0 ) {
      const <- -min.value + 10e-9 # then, min.value must be 10e-9
   } else {
      const <- 0 # no need to add any value 
   }
   pheno.vec.posit <- pheno.vec + const
   
   # box-cox trans
   if (lambda == 0) {
      pheno.trans <- log(pheno.vec.posit)
   } else {
      pheno.trans <- pheno.vec.posit ^ lambda
   }
   
   # substitue
   const.vec[trait] <- const
   df.toco.pheno.trans[[trait]] <- pheno.trans
}

# data frame of lambda and const
df.param <- data.frame("trait" = trait.all, "lambda" = lambda.vec, "const" = const.vec)

# save important objects
write.csv(df.toco.pheno, file = out.pheno, row.names = F)
write.csv(df.toco.pheno.trans, file = out.pheno.trans, row.names = F)
write.csv(df.info, file = out.info, row.names = F)
write.csv(df.param, file = out.BoxCoxParam, row.names = F)
