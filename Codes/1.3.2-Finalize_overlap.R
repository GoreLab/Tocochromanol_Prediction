# package
library(readxl)
library(ggplot2)

# load data
df.ID.match <- as.data.frame(read_xlsx("RESULT/1.3-CalcIbs/df.282.and.ames.ID.match_manual_curate.xlsx"))
IBS.mat <- as.matrix(read.csv("RESULT/1.3-CalcIbs/IbsMat.csv", row.names = 1))
AmesPheno <- read.csv("RESULT/1.1-MakeAmesPhenoData/AmesPheno.csv")
GbPheno <- read.csv("RESULT/1.2-MakeGbPhenoData/GbPheno.csv")

# Get 170 lines flagged as 'overlapped'
df.overlap.all <- df.ID.match[df.ID.match$In.Ames.manual == "yes", ]

# change text (e.g. AMES -> Ames) to match Ames IDs
df.overlap.all$GRIN.ID.no.space <- gsub("AMES", "Ames", df.overlap.all$GRIN.ID.no.space)
df.overlap.all$GRIN.ID.no.space <- gsub("CIZE", "CIze", df.overlap.all$GRIN.ID.no.space)

# Fix the GRIN ID for the three manually flagged lines
df.overlap.all$GRIN.ID.no.space[df.overlap.all$In.Ames == "no"] <- df.overlap.all$Ames.ID.manual[df.overlap.all$In.Ames == "no"]

# check
all(df.overlap.all$GBS.Sample.ID %in% rownames(IBS.mat)) # ok!
all(df.overlap.all$GRIN.ID.no.space %in% rownames(IBS.mat)) # ok!

# get the IBS between the same genotype in two panels
df.pair.IBS <- data.frame("Name.Matched" = df.overlap.all$Name.Matched, 
						  "ID.282" = df.overlap.all$GBS.Sample.ID,
						  "ID.Ames" = df.overlap.all$GRIN.ID.no.space,
						  "IBS" = NA)
for ( i in 1:nrow(df.pair.IBS) ) {
	id.row <- df.pair.IBS$ID.282[i]
	id.col <- df.pair.IBS$ID.Ames[i]
	ibs <- IBS.mat[id.row, id.col]
	df.pair.IBS$IBS[i] <- ibs
}
write.csv(df.pair.IBS, file = "RESULT/1.3-CalcIbs/df.self.IBS.for.170.lines.csv", row.names = F)

# show result
p <- ggplot(df.pair.IBS, aes(x = IBS))
p <- p + geom_histogram()
p <- p + ggtitle("IBS between the different GBS sample for the 'same' genotype")
ggsave(filename = "RESULT/1.3-CalcIbs/IBS_of_same_genotype_from_different_soruse.png", p,
	   width = 8, height = 6)

# identify the three low IBS lines & re-flag them: we treat them as different lines in our study
df.pair.IBS[df.pair.IBS$IBS < 0.8, ] # these four

# final set of the overlapped ones
df.overlap <- df.pair.IBS[df.pair.IBS$IBS >= 0.8, c("ID.282", "ID.Ames", "IBS")]


# make a data frame to summarize the line IDs used in this study
df.all.lines <- data.frame("Panel" = c(rep("Ames", nrow(AmesPheno)), rep("282", nrow(GbPheno))), 
						   "ID" = c(AmesPheno$ID, GbPheno$ID),
						   "Overlap" = NA,
						   "other.ID" = NA)
tf.dup <- df.all.lines$ID %in% c(df.overlap$ID.282, df.overlap$ID.Ames)
df.all.lines$Overlap[tf.dup] <- "Yes"
df.all.lines$Overlap[!tf.dup] <- "No"
m1 <- match(df.overlap$ID.282, df.all.lines$ID)
df.all.lines$other.ID[m1] <- df.overlap$ID.Ames
m2 <- match(df.overlap$ID.Ames, df.all.lines$ID)
df.all.lines$other.ID[m2] <- df.overlap$ID.282
write.csv(df.all.lines, "RESULT/1.3-CalcIbs/Summary_of_1702_lines.csv", row.names = F)











