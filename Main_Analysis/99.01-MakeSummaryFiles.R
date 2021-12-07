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
# ----- Baseline, AcrossPop
# ---------------------------------------------------------------------------- #
model.all <- c("BRR", "BayesB", "GBLUP", "GBLUP_12M")
train.all <- c("ames", "gb")
df.cor <- expand.grid("Model" = model.all, "Trait" = trait.all, "Train" = train.all, stringsAsFactors = F)
pred.mat <- NULL
for ( i in 1:nrow(df.cor) ) {
   # i-th case
   model <- df.cor$Model[i]
   trait <- df.cor$Trait[i]
   train <- df.cor$Train[i]
   
   # load
   if ( model %in% c("BRR", "BayesB") ) {
      f <- paste0("RESULT/2.2-PredictionAccuracy_trans/",
                  "yPred_trait=", trait, "_train=", train, "_model=", model, ".csv")
   }
   if ( model == "GBLUP" ) {
      f <- paste0("RESULT/2.3-GenPreFit_SinglePop_trans_GBLUP/",
                  "yPred_trait=", trait, "_train=", train, "_model=", model, ".csv")
   }
   if ( model == "GBLUP_12M" ) {
      f <- paste0("RESULT/2.4-GenPreFit_SinglePop_trans_GBLUP_12M/",
                  "yPred_trait=", trait, "_train=", train, "_model=GBLUP.csv")
   }
   pred.res <- fread(f, data.table = F)
   
   # check
   flag <- setequal(pred.res$ID, c(AmesPheno$ID, GbPheno$ID))
   if ( flag == FALSE ) { print("ERROR!!"); break }
   
   # observed values
   if ( train == "ames" ) { obs <- GbPheno[, trait]; names(obs) <- GbPheno$ID }
   if ( train == "gb" ) { obs <- AmesPheno[, trait]; names(obs) <- AmesPheno$ID }
   
   # predicted values (match data)
   m <- match(names(obs), pred.res$ID)
   pred <- pred.res$pred[m]
   pred.mat <- rbind(pred.mat, pred.res$pred)
   colnames(pred.mat) <- pred.res$ID
}
df.pred.all <- cbind(df.cor[, c("Train", "Trait", "Model")], pred.mat)
df.pred.all$Model[df.pred.all$Model == "GBLUP"] <- "GBLUP_341K"
df.pred.all$Model[df.pred.all$Model == "BRR"] <- "BRR_341K"
df.pred.all$Model[df.pred.all$Model == "BayesB"] <- "BayesB_341K"
df.pred.from.ames.to.282 <- df.pred.all[df.pred.all$Train == "ames", ]
df.pred.from.282.to.ames <- df.pred.all[df.pred.all$Train == "gb", ]
df.pred.from.ames.to.282.save <- df.pred.from.ames.to.282[, c("Trait", "Model", GbPheno$ID)]
df.pred.from.282.to.ames.save <- df.pred.from.282.to.ames[, c("Trait", "Model", AmesPheno$ID)]
write.csv(df.pred.from.ames.to.282.save, file = paste0(dir.save, "/BaselineModels_AmesTo282_BackTransPredVal.csv"), row.names = F)
write.csv(df.pred.from.282.to.ames.save, file = paste0(dir.save, "/BaselineModels_282ToAmes_BackTransPredVal.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- Baseline, WithinPop
# ---------------------------------------------------------------------------- #
model.all <- c("GBLUP", "BRR", "BayesB", "GBLUP_12M")
train.all <- c("ames", "gb")
df.cor <- expand.grid("Model" = model.all,
                      "Trait" = trait.all,
                      "Train" = train.all,
                      stringsAsFactors = F)
df.pred.Ames <- df.pred.282 <- df.all <- NULL
for ( i in 1:nrow(df.cor) ) {
   # i-th case
   model <- df.cor$Model[i]
   trait <- df.cor$Trait[i]
   train <- df.cor$Train[i]
   
   # prediction
   if ( model %in% c("BRR", "BayesB") ) {
      f <- paste0("RESULT/3.1-GenPreCv_SinglePop_trans/",
                  "yPred_Trait=", trait, "_Train=", train, "_Model=", model, ".csv")
   } 
   if ( model == "GBLUP" ) {
      f <- paste0("RESULT/3.2-GenPreCv_SinglePop_trans_GBLUP/",
                  "yPred_Trait=", trait, "_Train=", train, "_Model=", model, ".csv")
   }
   if ( model == "GBLUP_12M" ) {
      f <- paste0("RESULT/3.3-GenPreCv_SinglePop_trans_GBLUP_12M/",
                  "yPred_trait=", trait, "_train=", train, "_model=GBLUP.csv")
   }
   pred.res <- fread(f, data.table = F)
   
   # box-cox parameter
   if ( train == "ames" ) {
      lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
      const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   }
   if ( train == "gb" ) {
      lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
      const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
   }
   
   # observed values
   if ( train == "ames" ) { obs <- AmesPheno[, trait]; names(obs) <- AmesPheno$ID }
   if ( train == "gb" ) { obs <- GbPheno[, trait]; names(obs) <- GbPheno$ID }
   
   # predicted values (match data)
   m <- match(names(obs), pred.res$ID)
   pred.df <- pred.res[m, paste0("rep", formatC(1:10, width = 2, flag = "0"))]
   
   # accuracy
   reps <- colnames(pred.df)
   r.all <- c()
   pred.mat <- NULL
   for ( k in 1:length(reps) ) {
      name.k <- reps[k]
      pred <- pred.df[, name.k]
      
      # convert
      tmp.vec <- as.numeric(pred)
      if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
      tmp.vec.original <- tmp.vec.original - const
      pred.new <- tmp.vec.original
      names(pred.new) <- names(obs)
      
      # correlation
      r <- cor(obs, pred.new, use = "p")
      r.all <- c(r.all, r)
      
      # predicted values
      pred.mat <- rbind(pred.mat, pred.new)
   }
   
   # save predicted values
   df.pred.i <- data.frame("Model" = model, "Trait" = trait, "Train" = train,
                           "rep" = reps, pred.mat, row.names = NULL)
   if ( train == "ames" ) { df.pred.Ames <- rbind(df.pred.Ames, df.pred.i) }
   if ( train == "gb" ) { df.pred.282 <- rbind(df.pred.282, df.pred.i) }
}
colnames(df.pred.282)[5:ncol(df.pred.282)] <-  GbPheno$ID
df.pred.Ames$Train <- "Ames panel"
df.pred.282$Train <- "The 282 panel"
df.pred.Ames <- df.pred.Ames[, c(2, 1, 4, 5:ncol(df.pred.Ames))]
df.pred.282 <- df.pred.282[, c(2, 1, 4, 5:ncol(df.pred.282))]
df.pred.Ames$Model[df.pred.Ames$Model == "GBLUP"] <- "GBLUP_341K"; df.pred.282$Model[df.pred.282$Model == "GBLUP"] <- "GBLUP_341K"
df.pred.Ames$Model[df.pred.Ames$Model == "BRR"] <- "BRR_341K"; df.pred.282$Model[df.pred.282$Model == "BRR"] <- "BRR_341K"
df.pred.Ames$Model[df.pred.Ames$Model == "BayesB"] <- "BayesB_341K"; df.pred.282$Model[df.pred.282$Model == "BayesB"] <- "BayesB_341K"
colnames(df.pred.Ames)[colnames(df.pred.Ames) == "rep"] <- "Rep"
colnames(df.pred.282)[colnames(df.pred.282) == "rep"] <- "Rep"
write.csv(df.pred.Ames, file = paste0(dir.save, "/BaselineModels_WithinAmes_BackTransPredVal.csv"), row.names = F)
write.csv(df.pred.282, file = paste0(dir.save, "/BaselineModels_Within282_BackTransPredVal.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- NAM-QTL, AcrossPop
# ---------------------------------------------------------------------------- #
model.all <- c("250K", "1M", "SI")
train.all <- c("ames", "gb")
df.cor <- expand.grid("Model" = model.all, "Trait" = trait.all, "Train" = train.all, stringsAsFactors = F)
df.all <- pred.mat <- NULL
for ( i in 1:nrow(df.cor) ) {
   # i-th case
   model <- df.cor$Model[i]
   trait <- df.cor$Trait[i]
   train <- df.cor$Train[i]
   
   # load
   if ( train == "ames" ) {
      f <- paste0("RESULT/5.1-MultiKernel_Across_trans/",
                  "PredRes_AmesToGb_MultiKernel_", model, "_", trait, "_trans.csv")
      lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
      const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   }
   if ( train == "gb" ) {
      f <- paste0("RESULT/5.1-MultiKernel_Across_trans/",
                  "PredRes_GbToAmes_MultiKernel_", model, "_", trait, "_trans.csv")
      lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
      const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
   }
   pred.res <- fread(f, data.table = F)
   
   # check
   flag <- setequal(pred.res[, 1], c(AmesPheno$ID, GbPheno$ID))
   if ( flag == FALSE ) { print("ERROR!!"); break }
   
   # convert
   pred <- pred.res$UseAll
   tmp.vec <- as.numeric(pred)
   if (lambda == 0) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   pred.new <- tmp.vec.original
   names(pred.new) <- pred.res$GBS.ID
   
   # observed values
   if ( train == "ames" ) { obs <- GbPheno[, trait]; names(obs) <- GbPheno$ID }
   if ( train == "gb" ) { obs <- AmesPheno[, trait]; names(obs) <- AmesPheno$ID }
   
   # predicted values (match data)
   m <- match(names(obs), names(pred.new))
   pred.new.m <- pred.new[m]
   pred.mat <- rbind(pred.mat, pred.new)
}
df.pred.all <- cbind(df.cor[, c("Train", "Trait", "Model")], pred.mat)
df.pred.from.ames.to.282 <- df.pred.all[df.pred.all$Train == "ames", ]
df.pred.from.282.to.ames <- df.pred.all[df.pred.all$Train == "gb", ]
df.pred.from.ames.to.282.save <- df.pred.from.ames.to.282[, c("Trait", "Model", GbPheno$ID)]
df.pred.from.282.to.ames.save <- df.pred.from.282.to.ames[, c("Trait", "Model", AmesPheno$ID)]
write.csv(df.pred.from.ames.to.282.save, file = paste0(dir.save, "/UseNamQtl_AmesTo282_BackTransPredVal.csv"), row.names = F)
write.csv(df.pred.from.282.to.ames.save, file = paste0(dir.save, "/UseNamQtl_282ToAmes_BackTransPredVal.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- NAM-QTL, WithinPop
# ---------------------------------------------------------------------------- #
model.all <- c("250K", "1M", "SI")
pop.all <- c("ames", "gb")
df.summary <- expand.grid("Trait" = trait.all, "Model" = model.all, "Train" = pop.all, stringsAsFactors = F)
pred.df.all.ames <- pred.df.all.282 <- NULL
for ( i in 1:nrow(df.summary) ) {
   # parameter
   trait <- df.summary$Trait[i]
   model <- df.summary$Model[i]
   pop <- df.summary$Train[i]
   
   # correlation
   f <- paste0("RESULT/5.2-MultiKernel_Cv_trans/pred_", pop, "_", trait, "_", model, ".csv")
   pred.res <- fread(f, data.table = F)

   # save all predicted values
   pred.mat <- t(pred.res[, 3:ncol(pred.res)])
   colnames(pred.mat) <- pred.res$ID
   pred.df <- data.frame("Trait" = trait,
                         "Model" = model,
                         "Train" = pop,
                         "Rep" = rownames(pred.mat),
                         pred.mat, 
                         row.names = NULL)
   
   # stack the data frame of the predicted values
   if ( pop == "ames" ) { pred.df.all.ames <- rbind(pred.df.all.ames, pred.df) }
   if ( pop == "gb" ) { pred.df.all.282 <- rbind(pred.df.all.282, pred.df) }
}
colnames(pred.df.all.282)[5:ncol(pred.df.all.282)] <- GbPheno$ID
df.ames.save <- pred.df.all.ames[ c(1, 2, 4, 5:ncol(pred.df.all.ames))]
df.282.save <- pred.df.all.282[ c(1, 2, 4, 5:ncol(pred.df.all.282))]
write.csv(df.ames.save, file = paste0(dir.save, "/UseNamQtl_WithinAmes_BackTransPredVal.csv"), row.names = F)
write.csv(df.282.save, file = paste0(dir.save, "/UseNamQtl_Within282_BackTransPredVal.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- Transcriptome-based (witin-Ames only!)
# ---------------------------------------------------------------------------- #
# first five models
model.all <- c("ExpBLUP", "ExpBLUP_UseCandGene", 
               "GBLUP", "GBLUP+ExpBLUP", "GBLUP+ExpBLUP_UseCandGene")
df <- expand.grid("Model" = model.all, "Trait" = trait.all, stringsAsFactors = F)
df.pred <- NULL
for ( i in 1:nrow(df) ) {
   trait <- df$Trait[i]
   model <- df$Model[i]

   # prediction result
   f <- paste0("RESULT/4.1-GenExpPred_trans/pred_", trait, "_", model, ".csv")
   pred <- read.csv(f)
   
   # predicted values
   pred.mat <- t(pred[, -1])
   colnames(pred.mat) <- pred$ID
   pred.mat <- pred.mat[-1, ]
   df.pred.i <- data.frame("Trait" = trait,
                           "Model" = model,
                           "Rep" =  rownames(pred.mat),
                           pred.mat, row.names = NULL)
   df.pred <- rbind(df.pred, df.pred.i)
}
df.pred$Model[df.pred$Model == "ExpBLUP"] <- "TBLUP"
df.pred$Model[df.pred$Model == "ExpBLUP_UseCandGene"] <- "TBLUP_sub"
df.pred$Model[df.pred$Model == "GBLUP+ExpBLUP"] <- "GBLUP+TBLUP"
df.pred$Model[df.pred$Model == "GBLUP+ExpBLUP_UseCandGene"] <- "GBLUP+TBLUP_sub"
df.pred <- cbind.data.frame(df.pred[, c("Trait", "Model", "Rep")], 
                            df.pred[, !(colnames(df.pred) %in% c("Trait", "Model", "Rep"))])
# two additionl models
method.all <- c("MultiNamLargeGenesCovar", "MultiNamLargeGenes")
# loop for all traits
df.02 <- NULL
for ( tr in 1:length(trait.all) ) {
   # trait
   trait <- trait.all[tr]
   
   # methods
   for ( i in 1:length(method.all) ) {
      # method
      method <- method.all[i]
      
      # prediction result
      f <- paste0("MultiTraitModel/RESULT/CvRes_", trait, "_", method, ".csv")
      pred.df <- read.csv(f, row.names = 1)
      
      # back transform
      lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
      const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
      for ( j in 1:ncol(pred.df) ) {
         pred <- pred.df[, j]
         if (lambda == 0) { pred.ut.tmp <- exp(pred) } else { pred.ut.tmp <- pred ^ (1/lambda) }
         pred.ut <- pred.ut.tmp - const
         pred.df[, j] <- pred.ut
      }
      
      # predicted values
      pred.mat <- t(pred.df)
      colnames(pred.mat) <- rownames(pred.df)
      df.pred.i <- data.frame("Trait" = trait,
                              "Model" = method,
                              "Rep" =  rownames(pred.mat),
                              pred.mat, row.names = NULL)
      df.02 <- rbind(df.02, df.pred.i)
   }
}
flag <- all(colnames(df.pred) == colnames(df.02))
if ( flag ) { 
   df.pred.all <- rbind.data.frame(df.pred, df.02)
   write.csv(df.pred.all, file = paste0(dir.save, "/TrsctPred_WithinAmes_BackTransPredVal.csv"), row.names = F)
}

# ---------------------------------------------------------------------------- #
# ----- Genome features, AcrossPop
# ---------------------------------------------------------------------------- #
feature.all <- c("GERP", "MAF", "Prox_1kb", "Prox_10Kb", "Prox_100Kb", "RecRate")
train.all <- c("ames", "gb")
df.cor <- expand.grid("Feature" = feature.all, "Trait" = trait.all, "Train" = train.all, stringsAsFactors = F)
pred.all <- NULL
for ( i in 1:nrow(df.cor) ) {
   # i-th case
   feature <- df.cor$Feature[i]
   trait <- df.cor$Trait[i]
   train <- df.cor$Train[i]
   
   # prediction
   if ( train == "ames" ) {
      f <- paste0("RESULT/6.3-UseFeatures_AcrossPopPred_trans_12M/",
                  "PredRes_AmesToGb_MultiKernel_", feature, "_", trait, "_trans.csv")
   } else {
      f <- paste0("RESULT/6.3-UseFeatures_AcrossPopPred_trans_12M/",
                  "PredRes_GbToAmes_MultiKernel_", feature, "_", trait, "_trans.csv")
   }
   pred.res <- fread(f, data.table = F)
   
   # box-cox
   if ( train == "ames" ) {
      lambda <- BoxCoxParam.ames$lambda[BoxCoxParam.ames$trait == trait]
      const <- BoxCoxParam.ames$const[BoxCoxParam.ames$trait == trait]
   }
   if ( train == "gb" ) {
      lambda <- BoxCoxParam.gb$lambda[BoxCoxParam.gb$trait == trait]
      const <- BoxCoxParam.gb$const[BoxCoxParam.gb$trait == trait]
   }
   
   # observed values
   if ( train == "ames" ) { obs <- GbPheno[, trait]; names(obs) <- GbPheno$ID }
   if ( train == "gb" ) { obs <- AmesPheno[, trait]; names(obs) <- AmesPheno$ID }
   
   # convert
   pred <- pred.res$UseAll
   tmp.vec <- as.numeric(pred)
   if ( lambda == 0 ) { tmp.vec.original <- exp(tmp.vec) } else { tmp.vec.original <- tmp.vec ^ (1/lambda) }
   tmp.vec.original <- tmp.vec.original - const
   pred.new <- tmp.vec.original
   names(pred.new) <- pred.res$GBS.ID
   
   # match data
   m <- match(names(obs), names(pred.new))
   pred.new.m <- pred.new[m]
   
   # stack the predicted values
   pred.all <- rbind(pred.all, pred.new)
}
pred.all.df <- cbind(df.cor, pred.all)
pred.all.df$Feature[pred.all.df$Feature == "Prox_1kb"] <- "Prox_1Kb"
df.pred.from.ames.to.282 <- pred.all.df[pred.all.df$Train == "ames", ]
df.pred.from.282.to.ames <- pred.all.df[pred.all.df$Train == "gb", ]
df.pred.from.ames.to.282.save <- df.pred.from.ames.to.282[, c("Trait", "Feature", GbPheno$ID)]
df.pred.from.282.to.ames.save <- df.pred.from.282.to.ames[, c("Trait", "Feature", AmesPheno$ID)]
colnames(df.pred.from.ames.to.282.save)[colnames(df.pred.from.ames.to.282.save) == "Feature"] <- "Model"
colnames(df.pred.from.282.to.ames.save)[colnames(df.pred.from.282.to.ames.save) == "Feature"] <- "Model"
write.csv(df.pred.from.ames.to.282.save, file = paste0(dir.save, "/UseFeat_12M_AmesTo282_BackTransPredVal.csv"), row.names = F)
write.csv(df.pred.from.282.to.ames.save, file = paste0(dir.save, "/UseFeat_12M_282ToAmes_BackTransPredVal.csv"), row.names = F)

# ---------------------------------------------------------------------------- #
# ----- Genome features, WithinPop
# ---------------------------------------------------------------------------- #
model.all <- c("GERP", "MAF", "RecRate", "Prox_1kb", "Prox_10kb", "Prox_100kb")
pop.all <- c("ames", "gb")
df.summary <- expand.grid("Trait" = trait.all, "Model" = model.all, "Train" = pop.all,
                          stringsAsFactors = F)
pred.df.all.ames <- pred.df.all.282 <- NULL
for ( i in 1:nrow(df.summary) ) {
   # parameter
   trait <- df.summary$Trait[i]
   model <- df.summary$Model[i]
   pop <- df.summary$Train[i]
   
   # correlation
   f <- paste0("RESULT/6.4-UseFeatures_cv_12M/pred_", trait, "_", pop, "_GBLUP_Use_", model, "_cv.csv")
   pred.res <- fread(f, data.table = F)
   
   # save all predicted values
   pred.mat <- t(pred.res[, 3:ncol(pred.res)])
   colnames(pred.mat) <- pred.res$ID
   pred.df <- data.frame("Trait" = trait,
                         "Model" = model,
                         "Train" = pop,
                         "Rep" = rownames(pred.mat),
                         pred.mat, row.names = NULL)
   
   # stack the data frame of the predicted values
   if ( pop == "ames" ) { pred.df.all.ames <- rbind(pred.df.all.ames, pred.df) }
   if ( pop == "gb" ) { pred.df.all.282 <- rbind(pred.df.all.282, pred.df) }
}
colnames(pred.df.all.282)[5:ncol(pred.df.all.282)] <- GbPheno$ID
pred.df.all.ames$Model[pred.df.all.ames$Model == "Prox_1kb"] <- "Prox_1Kb"
pred.df.all.ames$Model[pred.df.all.ames$Model == "Prox_10kb"] <- "Prox_10Kb"
pred.df.all.ames$Model[pred.df.all.ames$Model == "Prox_100kb"] <- "Prox_100Kb"
pred.df.all.282$Model[pred.df.all.282$Model == "Prox_1kb"] <- "Prox_1Kb"
pred.df.all.282$Model[pred.df.all.282$Model == "Prox_10kb"] <- "Prox_10Kb"
pred.df.all.282$Model[pred.df.all.282$Model == "Prox_100kb"] <- "Prox_100Kb"
df.ames.save <- pred.df.all.ames[ c(1, 2, 4, 5:ncol(pred.df.all.ames))]
df.282.save <- pred.df.all.282[ c(1, 2, 4, 5:ncol(pred.df.all.282))]
write.csv(df.ames.save, file = paste0(dir.save, "/UseFeat_12M_WithinAmes_BackTransPredVal.csv"), row.names = F)
write.csv(df.282.save, file = paste0(dir.save, "/UseFeat_12M_Within282_BackTransPredVal.csv"), row.names = F)

