# package
library(data.table)
library(plotrix)

# mkdir
dir.save <- "RESULT/99-SummaryFiles"
dir.create(dir.save, recursive = T)

# load phenotype data
f.ames <- "RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv" # Ames Phenotype data
f.gb <- "RESULT/1.2-MakeGbPhenoData/GbPheno.csv" # Gb Phenotype data
AmesPheno <- fread(f.ames, data.table = F)
GbPheno <- fread(f.gb, data.table = F)

# load box-cox parameter
BoxCoxParam.ames <- read.csv("RESULT/1.1-MakeAmesPhenoData/Ames_BoxCoxParam.csv")
BoxCoxParam.gb <- read.csv("RESULT/1.2-MakeGbPhenoData/Gb_BoxCoxParam.csv")

# traits
trait.all <- c("a.T", "d.T", "g.T", "a.T3", "d.T3", "g.T3",
               "Total.Tocopherols", "Total.Tocotrienols", "Total.Tocochromanols")

# ---------------------------------------------------------------------------- #
# ----- NAM-QTL, AcrossPop: evalaution of each kernel
# ---------------------------------------------------------------------------- #
myFun.ChangeName <- function(x){
   num.lb <- regexpr("\\(", x)
   if ( num.lb == -1 ) {
      name.new <- x
   } else {
      name.new <- substr(x, 1, num.lb[1] - 1)
   }
   return(name.new)
}
MyFun.InvBoxCox <- function(x, lambda, const){
   if (lambda == 0) { x.orig <- exp(x) } else { x.orig <- x ^ (1/lambda) }
   x.orig.2 <- x.orig - const
   return(x.orig.2)
}
model.all <- c("250K", "1M", "SI")
scenario.all <- c("AmesToGb", "GbToAmes")
df.cor <- expand.grid("Model" = model.all, "Trait" = trait.all, "Scenario" = scenario.all, stringsAsFactors = F)
for ( i in 1:nrow(df.cor) ) {
   # parameter
   trait <- df.cor$Trait[i]
   scenario <- df.cor$Scenario[i]
   win <- df.cor$Model[i]
   if ( win == "250K" ) { win.v2 <- "250Kb" } else { win.v2 <- win }
   
   # Load Result 1. Fit MK-GBLUP
   f1 <- paste0("RESULT/5.1-MultiKernel_Across_trans/",
                "PredRes_", scenario, "_MultiKernel_", win, "_", trait, "_trans.csv")
   res1 <- fread(f1, data.table = F)
   res1 <- res1[, !(colnames(res1) %in% c("obs"))] # remove some columns (not needed)
   colnames(res1) <- sapply(colnames(res1), myFun.ChangeName, USE.NAMES = F)
   res1.mat <- as.matrix(res1[, -1])
   rownames(res1.mat) <- res1$GBS.ID
   
   # Load Result 2. Fit GBLUP (single-kernel) one-by-one for each kernel
   f2 <- paste0("RESULT/8.1-MultiKernel_Across_trans_eval_each_kern/",
                "PredRes_", scenario, "_EachKernel_", win.v2, "_", trait, "_trans.csv")
   res2 <- fread(f2, data.table = F)
   res2.mat <- as.matrix(res2[, -1])
   rownames(res2.mat) <- res2$GBS.ID
   
   # BC parameters
   if ( scenario == "AmesToGb" ) {
      L <- BoxCoxParam.ames[BoxCoxParam.ames$trait == trait, "lambda"]
      C <- BoxCoxParam.ames[BoxCoxParam.ames$trait == trait, "const"]
   }
   if ( scenario == "GbToAmes" ) {
      L <- BoxCoxParam.gb[BoxCoxParam.gb$trait == trait, "lambda"]
      C <- BoxCoxParam.gb[BoxCoxParam.gb$trait == trait, "const"]
   }
   
   # Inverse BC
   res1.mat.bt <- apply(res1.mat, 2, MyFun.InvBoxCox, lambda = L, const = C)
   res2.mat.bt <- apply(res2.mat, 2, MyFun.InvBoxCox, lambda = L, const = C)
   
   # BC predicted values
   df1 <- data.frame("Trait" = trait, 
                     "Model" = "MultiKernel",
                     "Window" = win, 
                     "Scenario" = scenario,
                     "QTL" = colnames(res1.mat.bt),
                     t(res1.mat.bt),
                     row.names = NULL)
   df2 <- data.frame("Trait" = trait, 
                     "Model" = "SingleKernel",
                     "Window" = win, 
                     "Scenario" = scenario,
                     "QTL" = colnames(res2.mat.bt),
                     t(res2.mat.bt),
                     row.names = NULL)
   
   # stack the data frame
   if ( i == 1) {
      df.MK <- df1
      df.SK <- df2      
   } else {
      df.MK <- rbind.data.frame(df.MK, df1)
      df.SK <- rbind.data.frame(df.SK, df2)
   }
   
   # name of accessions
   flag <- all(rownames(res1.mat.bt) == rownames(res2.mat.bt))   
   if ( flag == FALSE ) { print("ERROR!!"); break } else { acc.name <- rownames(res1.mat.bt) }
   if ( i == 1 ) {
      acc.name.init <- acc.name
   } else {
      flag.2 <- all(acc.name == acc.name.init)
      if ( flag == FALSE ) { print("ERROR!!"); break }
   }
}
df.NAM.PostHoc <- rbind.data.frame(df.SK, df.MK)
colnames(df.NAM.PostHoc)[6:ncol(df.NAM.PostHoc)] <- acc.name.init
df.NAM.PostHoc.FromAmesTo282 <- df.NAM.PostHoc[df.NAM.PostHoc$Scenario == "AmesToGb", ]
df.NAM.PostHoc.From282ToAmes <- df.NAM.PostHoc[df.NAM.PostHoc$Scenario == "GbToAmes", ]
df.NAM.PostHoc.FromAmesTo282.save <- df.NAM.PostHoc.FromAmesTo282[, c("Trait", "Model", "Window", "QTL", GbPheno$ID)]
df.NAM.PostHoc.From282ToAmes.save <- df.NAM.PostHoc.FromAmesTo282[, c("Trait", "Model", "Window", "QTL", AmesPheno$ID)]
write.csv(df.NAM.PostHoc.FromAmesTo282.save, file = paste0(dir.save, "/PostHoc_UseNamQtl_AmesTo282_EvalEachKernel.csv"), row.names = F)
write.csv(df.NAM.PostHoc.From282ToAmes.save, file = paste0(dir.save, "/PostHoc_UseNamQtl_282ToAmes_EvalEachKernel.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- Genome Feature, AcrossPop: evaluation of each kernel
# ---------------------------------------------------------------------------- #
model.all <- c("GERP", "MAF", "Prox_1kb", "Prox_10kb", "Prox_100kb", "RecRate")
scenario.all <- c("AmesToGb", "GbToAmes")
df.cor <- expand.grid("Model" = model.all,  "Trait" = trait.all, "Scenario" = scenario.all, stringsAsFactors = F)
for ( i in 1:nrow(df.cor) ) {
   # parameter
   model <- df.cor$Model[i]
   trait <- df.cor$Trait[i]
   scenario <- df.cor$Scenario[i]
   
   # Load Result 2. Fit GBLUP (single-kernel) one-by-one for each kernel
   f2 <- paste0("RESULT/8.5-FuncFeat_Across_trans_eval_each_kern_12M/",
                "PredRes_", scenario, "_EachKernel_", model, "_", trait, "_trans.csv")
   res2 <- fread(f2, data.table = F)
   res2.mat <- as.matrix(res2[, -1])
   rownames(res2.mat) <- res2$GBS.ID
   res2.mat <- res2.mat[, !(colnames(res2.mat) %in% "obs")]
   
   # Load Result 1. Fit MK-GBLUP
   f1 <- paste0("RESULT/6.3-UseFeatures_AcrossPopPred_trans_12M/",
                "PredRes_", scenario, "_MultiKernel_", model, "_", trait, "_trans.csv")
   res1 <- fread(f1, data.table = F)
   res1 <- res1[, !(colnames(res1) %in% c("obs"))] # remove some columns (not needed)
   colnames(res1) <- sapply(colnames(res1), myFun.ChangeName, USE.NAMES = F)
   res1.mat <- as.matrix(res1[, -1])
   rownames(res1.mat) <- res1$GBS.ID
   colnames(res1.mat)[-ncol(res1.mat)] <- colnames(res2.mat)
   
   # BC parameters
   if ( scenario == "AmesToGb" ) {
      L <- BoxCoxParam.ames[BoxCoxParam.ames$trait == trait, "lambda"]
      C <- BoxCoxParam.ames[BoxCoxParam.ames$trait == trait, "const"]
   }
   if ( scenario == "GbToAmes" ) {
      L <- BoxCoxParam.gb[BoxCoxParam.gb$trait == trait, "lambda"]
      C <- BoxCoxParam.gb[BoxCoxParam.gb$trait == trait, "const"]
   }
   
   # Inverse BC
   res1.mat.bt <- apply(res1.mat, 2, MyFun.InvBoxCox, lambda = L, const = C)
   res2.mat.bt <- apply(res2.mat, 2, MyFun.InvBoxCox, lambda = L, const = C)
   
   # BC predicted values
   df1 <- data.frame("Trait" = trait, 
                     "Model" = "MultiKernel",
                     "Feature" = model, 
                     "Scenario" = scenario,
                     "Grp" = colnames(res1.mat.bt),
                     t(res1.mat.bt),
                     row.names = NULL)
   df2 <- data.frame("Trait" = trait, 
                     "Model" = "SingleKernel",
                     "Feature" = model, 
                     "Scenario" = scenario,
                     "Grp" = colnames(res2.mat.bt),
                     t(res2.mat.bt),
                     row.names = NULL)
   
   # stack the data frame
   if ( i == 1 ) {
      df.MK <- df1
      df.SK <- df2      
   } else {
      flag <- all(colnames(df.SK) == colnames(df2))  
      if ( flag == FALSE ) { print("ERROR!!"); break }
      df.MK <- rbind.data.frame(df.MK, df1)
      df.SK <- rbind.data.frame(df.SK, df2)
   }
   
   # name of accessions
   flag <- all(rownames(res1.mat.bt) == rownames(res2.mat.bt))   
   if ( flag == FALSE ) { print("ERROR!!"); break } else { acc.name <- rownames(res1.mat.bt) }
   if ( i == 1 ) {
      acc.name.init <- acc.name
   } else {
      flag.2 <- all(acc.name == acc.name.init)
      if ( flag == FALSE ) { print("ERROR!!"); break }
   }
}
df.Feat.PostHoc <- rbind.data.frame(df.SK, df.MK)
colnames(df.Feat.PostHoc)[6:ncol(df.Feat.PostHoc)] <- acc.name.init
df.Feat.PostHoc$Feature[df.Feat.PostHoc$Feature == "Prox_1kb"] <- "Prox_1Kb"
df.Feat.PostHoc$Feature[df.Feat.PostHoc$Feature == "Prox_10kb"] <- "Prox_10Kb"
df.Feat.PostHoc$Feature[df.Feat.PostHoc$Feature == "Prox_100kb"] <- "Prox_100Kb"
colnames(df.Feat.PostHoc)[colnames(df.Feat.PostHoc) == "Grp"] <- "Group"
df.Feat.PostHoc.FromAmesTo282 <- df.Feat.PostHoc[df.Feat.PostHoc$Scenario == "AmesToGb", ]
df.Feat.PostHoc.From282ToAmes <- df.Feat.PostHoc[df.Feat.PostHoc$Scenario == "GbToAmes", ]
df.Feat.PostHoc.FromAmesTo282.save <- df.Feat.PostHoc.FromAmesTo282[, c("Trait", "Model", "Feature", "Group", GbPheno$ID)]
df.Feat.PostHoc.From282ToAmes.save <- df.Feat.PostHoc.From282ToAmes[, c("Trait", "Model", "Feature", "Group", AmesPheno$ID)]
write.csv(df.Feat.PostHoc.FromAmesTo282.save, file = paste0(dir.save, "/PostHoc_UseFeat_AmesTo282_EvalEachKernel_BackTransPredVal.csv"), row.names = F)
write.csv(df.Feat.PostHoc.From282ToAmes.save, file = paste0(dir.save, "/PostHoc_UseFeat_282ToAmes_EvalEachKernel_BackTransPredVal.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- Transcriptome-based prediction with different number of PEER factors
# ---------------------------------------------------------------------------- #
model.all <- c("ExpBLUP", "ExpBLUP_UseCandGene", 
               "GBLUP+ExpBLUP", "GBLUP+ExpBLUP_UseCandGene")
peer.all <- paste0("K", 1:25)
df <- expand.grid("Model" = model.all, "Trait" = trait.all, "PEER" = peer.all, stringsAsFactors = F)
df.pred <- NULL
for ( i in 1:nrow(df) ) {
   trait <- df$Trait[i]
   model <- df$Model[i]
   peer <- df$PEER[i]
   
   # prediction result
   f <- paste0("RESULT/4.3-GenExpPred_trans_TestPeer/pred_", trait, "_", model, "_UsePeer_", peer, ".csv")
   pred <- read.csv(f)
   
   # predicted values
   pred.mat <- t(pred[, -1])
   colnames(pred.mat) <- pred$ID
   pred.mat <- pred.mat[-1, ]
   df.pred.i <- data.frame("Trait" = trait,
                           "Model" = model,
                           "PEER" = peer,
                           "Rep" =  rownames(pred.mat),
                           pred.mat, row.names = NULL)
   df.pred <- rbind(df.pred, df.pred.i)
}
df.pred$Model[df.pred$Model == "ExpBLUP"] <- "TBLUP"
df.pred$Model[df.pred$Model == "ExpBLUP_UseCandGene"] <- "TBLUP_sub"
df.pred$Model[df.pred$Model == "GBLUP+ExpBLUP"] <- "GBLUP+TBLUP"
df.pred$Model[df.pred$Model == "GBLUP+ExpBLUP_UseCandGene"] <- "GBLUP+TBLUP_sub"
df.pred$PEER <- as.integer(gsub("K", "", df.pred$PEER))
colnames(df.pred)[colnames(df.pred) == "PEER"] <- "n.PEER.Factor"
write.csv(df.pred, file = paste0(dir.save, "/PostHoc_TrsctPred_WithinAmes_EvalNumPeerFact_BackTransPredVal.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- GBLUP with random X SNPs: From Ames to 282
# ---------------------------------------------------------------------------- #
snp.all <- as.character(as.integer(c(500000, 1000000, 5000000)))
dir.in <- "RESULT/9.2-GBLUP_with_random_SNPs_AcrossPop_12M/"
df <- expand.grid("Trait" = trait.all, "n.SNP" = snp.all, stringsAsFactors = F)
for ( i in 1:nrow(df) ) {
   # param
   trait <- df$Trait[i]
   n.snp <- df$n.SNP[i]
   
   # load file
   file.in <- paste0("PredRes_AmesToGb_RandomSnps_", n.snp, "_", trait, "_trans.csv")
   pred.res <- fread(file = paste0(dir.in, file.in), data.table = F)
   tf <- colnames(pred.res) %in% c("GBS.ID", paste0("Rep", formatC(1:10, width = 2 , flag = "0")))
   pred.res <- pred.res[, tf]
   
   # back-trans
   lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
   const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   for ( j in 2:ncol(pred.res) ) {
      pred <- pred.res[, j]
      tmp.vec <- as.numeric(pred)
      if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
      tmp.vec.original <- tmp.vec.original - const
      pred.new <- tmp.vec.original
      pred.res[, j] <- pred.new # substitute the converted values
   }
   
   # obs 
   obs <- GbPheno[, trait]; names(obs) <- GbPheno$ID
   
   # correlation
   m <- match(names(obs), pred.res$GBS.ID)
   pred.res.m <- pred.res[m, ]
   pred.res.m.t <-  t(pred.res.m[, -1])
   colnames(pred.res.m.t) <- pred.res.m[, 1]
   df.pred <- data.frame("Trait" = trait, "n.SNP" = n.snp, "Rep" = rownames(pred.res.m.t), pred.res.m.t)
   if ( i == 1 ) { df.pred.all <- df.pred } else { df.pred.all <- rbind.data.frame(df.pred.all, df.pred) }
}
rownames(df.pred.all) <- NULL
colnames(df.pred.all)[4:ncol(df.pred.all)] <- GbPheno$ID
write.csv(df.pred.all, file = paste0(dir.save, "/PostHoc_UseRandomSnpFrom12M_AmesTo282_BackTransPredVal.csv"), row.names = F)


# ---------------------------------------------------------------------------- #
# ----- GBLUP with random X SNPs: From 282 to Ames
# ---------------------------------------------------------------------------- #
snp.all <- as.character(as.integer(c(500000, 1000000, 5000000)))
df <- expand.grid("Trait" = trait.all, "n.SNP" = snp.all, stringsAsFactors = F)
df.accuracy <- NULL
for ( i in 1:nrow(df) ) {
   # param
   trait <- df$Trait[i]
   n.snp <- df$n.SNP[i]
   
   # load file
   file.in <- paste0("PredRes_GbToAmes_RandomSnps_", n.snp, "_", trait, "_trans.csv")
   pred.res <- fread(file = paste0(dir.in, file.in), data.table = F)
   tf <- colnames(pred.res) %in% c("GBS.ID", paste0("Rep", formatC(1:10, width = 2 , flag = "0")))
   pred.res <- pred.res[, tf]
   
   # back-trans
   lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
   const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
   for ( j in 2:ncol(pred.res) ) {
      pred <- pred.res[, j]
      tmp.vec <- as.numeric(pred)
      if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
      tmp.vec.original <- tmp.vec.original - const
      pred.new <- tmp.vec.original
      pred.res[, j] <- pred.new # substitute the converted values
   }
   
   # obs 
   obs <- AmesPheno[, trait]; names(obs) <- AmesPheno$ID
   
   # correlation
   m <- match(names(obs), pred.res$GBS.ID)
   pred.res.m <- pred.res[m, ]
   pred.res.m.t <-  t(pred.res.m[, -1])
   colnames(pred.res.m.t) <- pred.res.m[, 1]
   df.pred <- data.frame("Trait" = trait, "n.SNP" = n.snp, "Rep" = rownames(pred.res.m.t), pred.res.m.t)
   if ( i == 1 ) { df.pred.all <- df.pred } else { df.pred.all <- rbind.data.frame(df.pred.all, df.pred) }
}
rownames(df.pred.all) <- NULL
colnames(df.pred.all)[4:ncol(df.pred.all)] <- AmesPheno$ID
write.csv(df.pred.all, file = paste0(dir.save, "/PostHoc_UseRandomSnpFrom12M_282ToAmes_BackTransPredVal.csv"), row.names = F)





