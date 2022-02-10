# Make Gb phenotype Data

# source
library(gdata)
library(data.table)

# define object
trait.all <- c("aT", "dT", "gT",
               "aT3", "dT3", "gT3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")
trait.all.rename <- c("a.T", "d.T", "g.T",
                      "a.T3", "d.T3", "g.T3",
                      "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# file input
file.pheno.all <- "RAWDATA/Data_282/Lipka_2013_G3_TableS1_tocochromanol_BLUP.xlsx"
file.indv <- "RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv"
file.info <- "RAWDATA/Data_282/282.taxa.list.final.2020.02.06_full.csv"

# file output
dir.out <- "RESULT/1.2-MakeGbPhenoData"
dir.out.bc.log <- "RESULT/1.2-MakeGbPhenoData/BoxCoxLog"
out.pheno <- "RESULT/1.2-MakeGbPhenoData/GbPheno.csv"
out.pheno.trans <- "RESULT/1.2-MakeGbPhenoData/GbPheno_trans.csv"
out.info <- "RESULT/1.2-MakeGbPhenoData/GbInfo.csv"
out.BoxCoxParam <- "RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv"

# mkdir
dir.create(dir.out, recursive = TRUE)
dir.create(dir.out.bc.log, recursive = TRUE)

# load phenotype file
pheno.all <- read.xls(file.pheno.all, stringsAsFactors = F, skip = 1)
pheno.all <- pheno.all[, c("Sample.ID", trait.all)] # remove columns not needed
colnames(pheno.all) <- c("Sample.ID", trait.all.rename)

# load file of 242 accessions to be kept
info <- read.csv(file.info, stringsAsFactors = F)

# phenotype data to be kept
tf <- pheno.all$Sample.ID %in% info$Name.Lipka
table(tf) # TRUE should be 242 -- OK!!

# 242 phenotype datas
df.toco.pheno <- pheno.all[tf, ]

# load ID of the genotype data
indv.tmp <- fread(file.indv, header = FALSE)
id.geno <- indv.tmp$V1

# check if all phenotyped accessions are in the genotye data
setdiff(info$GBS.taxa, id.geno) # perfect!

# attach GBS name to the phenotype data frame
all.equal(df.toco.pheno$Sample.ID, info$Name.Lipka) # check -> OK
df.toco.pheno$Sample.ID <- info$GBS.taxa
colnames(df.toco.pheno)[1] <- "ID"

# ------------------------------------------------------------------------------------------------- #
# Use boxcox function to find the optimal lambda
source("0.1-MyFun_BoxcoxFunction.R")
lambda.vec <- rep(NA, times = length(trait.all.rename))
names(lambda.vec) <- trait.all.rename
for (tr in 1:length(trait.all.rename)) {
   trait <- trait.all.rename[tr]
   res.list <- boxcox.fitting(pheno = pheno.all, 
                              curr.trait = trait,
                              filename1 = paste0(dir.out.bc.log, "/", trait, "_lambda_values.jpeg"),
                              filename2 = paste0(dir.out.bc.log, "/", trait, "_Distrib_After_Box_plot.jpeg"))
   lambda.vec[trait] <- res.list$LAMBDA
}

# make transformed phenotype data
const.vec <- lambda.vec # set a object
df.toco.pheno.trans <- df.toco.pheno # set a object
for (trait in trait.all.rename) {
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
df.param <- data.frame("trait" = trait.all.rename, "lambda" = lambda.vec, "const" = const.vec)

# save important objects
write.csv(df.toco.pheno, file = out.pheno, row.names = F)
write.csv(df.toco.pheno.trans, file = out.pheno.trans, row.names = F)
write.csv(info, file = out.info, row.names = F)
write.csv(df.param, file = out.BoxCoxParam, row.names = F)
