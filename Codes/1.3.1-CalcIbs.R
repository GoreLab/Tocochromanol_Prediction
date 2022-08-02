# source
library(data.table)
library(readxl)

# mkdir 
dir.create("RESULT/1.3-CalcIbs")


# ---------------------------------------------------------------------------- #
# ----- IBS calculation
# ---------------------------------------------------------------------------- #
# load geotype & map
indv.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("RAWDATA/Genotype/merged_ames_282_LD0.1_all_chr.012.pos")
map <- data.frame("chr" = map.tmp$V1, "pos" = map.tmp$V2)
rm(indv.tmp); rm(geno.tmp); rm(map.tmp); gc(); gc() # remove non-essential object

# calculate pairwise IBS
P <- ncol(geno.mat)
score.mat <- geno.mat - 1
XtX <- tcrossprod(score.mat)
IBS.mat <- (XtX + P) / (2 * P)

# save the matrix
write.csv(IBS.mat, "RESULT/1.3-CalcIbs/IbsMat.csv")


# ---------------------------------------------------------------------------- #
# ----- Compare Ames and 282 and find overlap
# ---------------------------------------------------------------------------- #
# load phenotype data
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv")
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv")

# load ID match file
df.ID <- read_xlsx("RAWDATA/tpg220160-sup-0001-tables.xlsx", sheet = 1, skip = 1)
df.ID <- as.data.frame(df.ID)

# make new ID (remove GBS sample numbers)
x <- GbPheno$ID[1]
myfun.ConvertName <- function (x) {
	vec.x <- strsplit(x, ":")[[1]]
	new.x <- paste(vec.x[-length(vec.x)], sep = "_")
	return(new.x)
}
GbPheno$NewID <- sapply(GbPheno$ID, FUN = myfun.ConvertName, USE.NAMES = F)

# (1) Make a data frame for the accessions with perfect match
id.perfect.match <- intersect(GbPheno$NewID, df.ID$`Inbred line`)
m1 <- match(id.perfect.match, GbPheno$NewID)
m2 <- match(id.perfect.match, df.ID$`Inbred line`)
df.p.match <- data.frame("GBS.Sample.ID" = GbPheno$ID[m1],
						 "GBS.Sample.Name" = GbPheno$NewID[m1],
						 "Name.Matched" = df.ID$`Inbred line`[m2],
						 "Matching.Status" = "perfect match",
						 "GRIN.ID.no.space" = df.ID$`GRIN ID no space`[m2],
						 "Name.in.AllZeaGBSv2.7" = df.ID$`Name Found in AllZeaGBSv2.7 (Glaubitz et al. (2014) PLoS ONE 9:e90346)`[m2],
						 "Name.in.TPJ.study" = df.ID$`Name Used in this study`[m2],
						 "Number.of.GBS.Samples.in.TPJ.study" = df.ID$`Number of GBS Samples`[m2],
						 "Population.in.TPJ.study" = df.ID$Population[m2])

# (2) Make a data frame for the accessions using capitalized name
id.remain <- setdiff(GbPheno$NewID, id.perfect.match)
df.ID.sub <- df.ID[!(df.ID$`Inbred line` %in% id.perfect.match), ]
df.ID.sub.train <- df.ID.sub[df.ID.sub$Population == "Training", ]
m1 <- match(id.remain, GbPheno$NewID)
m2 <- match(toupper(id.remain), df.ID.sub.train$`Inbred line`)
df.u.match <- data.frame("GBS.Sample.ID" = GbPheno$ID[m1],
						 "GBS.Sample.Name" = GbPheno$NewID[m1],
						 "Name.Matched" = df.ID.sub.train$`Inbred line`[m2],
						 "Matching.Status" = "Uppercase match",
						 "GRIN.ID.no.space" = df.ID.sub.train$`GRIN ID no space`[m2],
						 "Name.in.AllZeaGBSv2.7" = df.ID.sub.train$`Name Found in AllZeaGBSv2.7 (Glaubitz et al. (2014) PLoS ONE 9:e90346)`[m2],
						 "Name.in.TPJ.study" = df.ID.sub.train$`Name Used in this study`[m2],
						 "Number.of.GBS.Samples.in.TPJ.study" = df.ID.sub.train$`Number of GBS Samples`[m2],
						 "Population.in.TPJ.study" = df.ID.sub.train$Population[m2])

# merge the two data frames
df.match <- rbind.data.frame(df.p.match, df.u.match)

# get GRIN ID
GRIN.ID.282 <- df.match$GRIN.ID.no.space
GRIN.ID.Ames <- toupper(AmesPheno$ID)
GRIN.ID.overlap <- intersect(GRIN.ID.282, GRIN.ID.Ames)
df.match$In.Ames <- "no"
df.match$In.Ames[df.match$GRIN.ID.no.space %in% GRIN.ID.overlap] <- "yes"

# max.IBS 
IBS.mat.across.pop <- IBS.mat[df.match$GBS.Sample.ID, AmesPheno$ID]
maxIBS <- apply(IBS.mat.across.pop, 1, max)
w.maxIBS <- apply(IBS.mat.across.pop, 1, which.max)
Ames.acc.w.maxIBS <- colnames(IBS.mat.across.pop)[w.maxIBS]
df.IBS <- data.frame("GBS.ID" = names(maxIBS),
					 "max.IBS" = maxIBS,
					 "Ames.Accession.with.max.IBS" = Ames.acc.w.maxIBS,
					 row.names = NULL)
m <- match(df.match$GBS.Sample.ID, df.IBS$GBS.ID)
df.IBS.all <- cbind.data.frame(df.match, 
							   "max.IBS" = df.IBS[m, "max.IBS"],
							   "Ames.Accession.with.max.IBS" = df.IBS[m, "Ames.Accession.with.max.IBS"])
df.IBS.sub <- df.IBS.all[, c("GBS.Sample.ID", "Name.Matched", 
							 "Matching.Status", "GRIN.ID.no.space", "Number.of.GBS.Samples.in.TPJ.study",
							 "In.Ames", "max.IBS", "Ames.Accession.with.max.IBS")]
write.csv(df.IBS.sub, file = "RESULT/1.3-CalcIbs/df.282.and.ames.ID.match.csv")

# 
# table(df.IBS.sub$In.Ames, df.IBS.sub$max.IBS < 0.95)
# summary(df.IBS$max.IBS[df.IBS.sub$In.Ames == "no"])
# 
# 
# IBS.within.Ames <- IBS.mat[AmesPheno$ID, AmesPheno$ID]
# IBS.within.282 <- IBS.mat[GbPheno$ID, GbPheno$ID]
# summary(IBS.within.Ames[upper.tri(IBS.within.Ames)])
# summary(IBS.within.282[upper.tri(IBS.within.282)])
# hist(IBS.within.Ames[upper.tri(IBS.within.Ames)])
# hist(IBS.within.282[upper.tri(IBS.within.282)])
# 
# diag(IBS.within.Ames) <- 0 # heuristic
# diag(IBS.within.282) <- 0  # heuristic
# maxIBS.within.Ames <- apply(IBS.within.Ames, 1, max)
# maxIBS.within.282 <- apply(IBS.within.282, 1, max)
# summary(maxIBS.within.Ames)
# summary(maxIBS.within.282)
# hist(maxIBS.within.Ames)
# hist(maxIBS.within.282)
# 
# 
# summary(maxIBS.within.Ames)
# quantile(maxIBS.within.Ames, 0.7)
# sum(maxIBS.within.Ames > 0.99) / length(maxIBS.within.Ames)
# 
# which(maxIBS.within.Ames > 0.99)
# 
# 
# which.max(IBS.within.Ames["Ames19316", ])
# 
# # #####
# # table(df.IBS.sub$In.Ames)
# # df.IBS.sub[df.IBS.sub$In.Ames == "no", ]
# # summary(df.IBS.sub$max.IBS[df.IBS.sub$In.Ames == "no"])
# # hist(df.IBS.sub$max.IBS[df.IBS.sub$In.Ames == "no"])
# # table(df.IBS.sub$max.IBS > 0.90 & df.IBS.sub$In.Ames == "no")
# # 
# # 
# # df.non.ames.tmp <- df.IBS.sub[df.IBS.sub$In.Ames == "no", ]
# # df.non.ames.tmp.ord <- df.non.ames.tmp[order(df.non.ames.tmp$max.IBS, decreasing = TRUE), ]
# # head(df.non.ames.tmp.ord)
# # #####
# # 
# # 
# # # change column names
# # df.282.sample <- data.frame("GBS.ID.282" = df.IBS.sub$GBS.Sample.ID,
# # 							"GRIN.ID.in.Dzievit.et.al" = df.IBS.sub$GRIN.ID.no.space,
# # 							"GRIN.ID.manual" = NA,
# # 							"Is.in.Ames" = NA,
# # 							"IBS" = NA,
# # 							"Overlap" = NA)
# # manual.id <- df.282.sample$GRIN.ID.in.Dzievit.et.al
# # for ( i in 1:length(manual.id) ) {
# # 	x <- manual.id[i]
# # 	x <- gsub("AMES", "Ames", x)
# # 	x <- gsub("CIZE", "CIze", x)
# # 	manual.id[i] <- x
# # }
# # df.282.sample$GRIN.ID.manual <- manual.id
# # df.282.sample$GRIN.ID.manual[df.282.sample$GBS.ID.282 == "H49:250040032"] <- "Ames26787"
# # df.282.sample$GRIN.ID.manual[df.282.sample$GBS.ID.282 == "CI64:250040128"] <- "CIze64"
# # df.282.sample$GRIN.ID.manual[df.282.sample$GBS.ID.282 == "Ab28A:250040159"] <- "Ames18999"
# # 
# # # In Ames?
# # tf <- df.282.sample$GRIN.ID.manual %in% colnames(IBS.mat.across.pop)
# # df.282.sample$Is.in.Ames[tf] <- "Yes"
# # df.282.sample$Is.in.Ames[!tf] <- "No"
# # 
# # # IBS
# # for ( i in 1:nrow(df.282.sample) ) {
# # 	if ( df.282.sample$Is.in.Ames[i] == "Yes" ) {
# # 		id.282 <- df.282.sample$GBS.ID.282[i]
# # 		id.ames <- df.282.sample$GRIN.ID.manual[i]
# # 		df.282.sample$IBS[i] <- IBS.mat.across.pop[id.282, id.ames]
# # 	} 
# # }
# # 
# # hist(df.282.sample$IBS)
# # sort(df.282.sample$IBS)
# # 
# # # overlap
# # tf <- (df.282.sample$Is.in.Ames == "Yes") & (df.282.sample$IBS > 0.8)
# # df.282.sample$Overlap[tf] <- "Overlap"
# # df.282.sample$Overlap[!tf] <- "Non-Overlap"
# # 
# # # 
# # png(filename = paste0(dir.save, "/Fig.self.IBS.hitogram.png"), width = 500, height = 400)
# # hist(df.282.sample$IBS, breaks = 20, main = "Histogram of self-IBS",
# # 	 xlab = "pairwise IBS")
# # dev.off()
# # 
# # # write table 
# # write.csv(df.282.sample, file = paste0(dir.save, "/Overlap_data.csv"), row.names = F)
# # 
# # 
# # ibs <- fread("C:/Users/rt475/Desktop/AP282_ibs.mibs")
# # ibs.id <- fread("C:/Users/rt475/Desktop/AP282_ibs.mibs.id", header = F)
# # ibs[1:5, 1:5]
# # ibs <- as.matrix(ibs)
# # rownames(ibs) <- rownames(IBS.mat)
# # rownames(ibs) <- colnames(IBS.mat)
# # 
# # head(ibs.id)
# # ibs.id[which(ibs.id$V1 != ibs.id$V2), ]
# # IBS.mat[which(ibs.id$V1 != ibs.id$V2), which(ibs.id$V1 != ibs.id$V2)]
# # 
# # 
# # A <- ibs
# # B <- round(IBS.mat, 6)
# # plot(A[upper.tri(A)], B[upper.tri(B)])
# # 
# # 
# # dalta <- A[upper.tri(A)] - B[upper.tri(B)]
# # hist(dalta)
# # max(dalta)

