# source
library(data.table)
library(readxl)

# mkdir 
dir.create("RESULT/1.3-CalcIbs", recursive = T)

# ---------------------------------------------------------------------------- #
# ----- IBS calculation
# ---------------------------------------------------------------------------- #
# load geotype & map
indv.tmp <- fread("merged_ames_282_LD0.1_all_chr.012.indv", header = FALSE)
geno.tmp <- fread("merged_ames_282_LD0.1_all_chr.012")
geno.mat <- as.matrix(geno.tmp[, -1]) # 1st col is ID.
rownames(geno.mat) <- indv.tmp$V1
map.tmp <- fread("merged_ames_282_LD0.1_all_chr.012.pos")
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
# ----- Check IBS
# ---------------------------------------------------------------------------- #
#
IBS.mat <- read.csv("RESULT/1.3-CalcIbs/IbsMat.csv")
IBS.mat <- as.matrix(IBS.mat[, -1])
rownames(IBS.mat) <- colnames(IBS.mat)

#
df.LineName <- data.table::fread("Table S1.txt", data.table = F)
colnames(df.LineName) <- c("Panel", "Inbred.line.name", "overlap")
name.ames <- df.LineName$Inbred.line.name[df.LineName$Panel == "Ames"]
name.good <- df.LineName$Inbred.line.name[df.LineName$Panel == "Goodman"]
overlapping.names <- intersect(name.ames, name.good)

#
df.ames.info <- data.frame(read_xlsx("vitamaize_accession_removal_20191113.xlsx"))
all(overlapping.names %in% df.ames.info$Pedigree_vita) # OK
ames.id.overlap <- df.ames.info$Accession.Number_ion[match(overlapping.names, df.ames.info$Pedigree_vita)]

# pairwise IBS
colnames(IBS.mat)[1463:1704] <- rownames(IBS.mat)[1463:1704] <- substr(rownames(IBS.mat)[1463:1704], 1, nchar(rownames(IBS.mat)[1463:1704]) - 10)

# change some names of Goodman to match with GBS samples
overlapping.names <- gsub(" ", "", overlapping.names)
overlapping.names <- gsub("-", ".", overlapping.names)
overlapping.names[overlapping.names == "38.11"] <- "X38.11"
overlapping.names[overlapping.names == "33.16"] <- "X33.16"
overlapping.names[overlapping.names == "4226"] <- "X4226"
overlapping.names[overlapping.names == "Oh7B"] <- "OH7B"
overlapping.names[overlapping.names == "Ms71"] <- "MS71"
overlapping.names[overlapping.names == "Va102"] <- "VA102"
overlapping.names[overlapping.names == "DE2"] <- "DE_2"
overlapping.names[overlapping.names == "DE3"] <- "DE_3"

# 
m.row <- match(ames.id.overlap, rownames(IBS.mat))
m.col <- match(overlapping.names, colnames(IBS.mat))
pair.ibs <- diag(IBS.mat[m.row, m.col])

#
pair.ibs[pair.ibs < 0.8]
overlapping.names[pair.ibs < 0.8] # four lines with IBS < 0.8

# # ---------------------------------------------------------------------------- #
# # ----- Another Check (considering the duplication within Ames)
# # ---------------------------------------------------------------------------- #
# # load IBS matrix
# IBS.mat <- read.csv("RESULT/1.3-CalcIbs/IbsMat.csv")
# IBS.mat <- as.matrix(IBS.mat[, -1])
# rownames(IBS.mat) <- colnames(IBS.mat)
# colnames(IBS.mat)[1463:1704] <- rownames(IBS.mat)[1463:1704] <- substr(rownames(IBS.mat)[1463:1704], 1, nchar(rownames(IBS.mat)[1463:1704]) - 10)
# 
# # load table S1 (with Ames Accession ID for the ones duplicated within Ames)
# df.LineName <- as.data.frame(read_excel("Table S1.xlsx"))
# 
# # Overlap between panels
# name.ames <- df.LineName$Inbred.line.name[df.LineName$Panel == "Ames"]
# name.good <- df.LineName$Inbred.line.name[df.LineName$Panel == "Goodman"]
# overlapping.names <- intersect(name.ames, name.good)
# ames.id.overlap <- df.ames.info$Accession.Number_ion[match(overlapping.names, df.ames.info$Pedigree_vita)]
# df.tmp <- data.frame("name" = overlapping.names,
#                      "id" = ames.id.overlap)
# 
# # change Goodman names to match with GBS samples
# overlapping.names.new <- gsub(" ", "", overlapping.names)
# overlapping.names.new <- gsub("-", ".", overlapping.names.new)
# overlapping.names.new[overlapping.names.new == "38.11"] <- "X38.11"
# overlapping.names.new[overlapping.names.new == "33.16"] <- "X33.16"
# overlapping.names.new[overlapping.names.new == "4226"] <- "X4226"
# overlapping.names.new[overlapping.names.new == "Oh7B"] <- "OH7B"
# overlapping.names.new[overlapping.names.new == "Ms71"] <- "MS71"
# overlapping.names.new[overlapping.names.new == "Va102"] <- "VA102"
# overlapping.names.new[overlapping.names.new == "DE2"] <- "DE_2"
# overlapping.names.new[overlapping.names.new == "DE3"] <- "DE_3"
# df.tmp$name.new <- overlapping.names.new
# 
# # Check IBS
# df.low.ibs <- NULL
# for ( i in 1:length(overlapping.names) ) {
#   name <- overlapping.names[i]
#   df.LineName.sub <- df.LineName[df.LineName$Inbred.line.name == name, ]
#   
#   if ( nrow(df.LineName.sub) == 3 ) {
#     id.ames <- df.LineName.sub$Accession.Name.or.Site.Identifier[!is.na(df.LineName.sub$Accession.Name.or.Site.Identifier)]
#     id.goodman <- df.tmp$name.new[df.tmp$name == name]
#     ibs <- IBS.mat[id.goodman, id.ames]
#     if ( max(ibs) < 0.8 ) { print(name); df.low.ibs <- rbind(df.low.ibs, data.frame("Ames" = names(ibs)[which.min(ibs)], "Goodman" = id.goodman)) }
#   } else {
#     id.ames <- df.tmp$id[df.tmp$name == name]
#     id.goodman <- df.tmp$name.new[df.tmp$name == name]
#     ibs <- IBS.mat[id.goodman, id.ames]
#     if ( ibs < 0.8 ) { print(name); df.low.ibs <- rbind(df.low.ibs, data.frame("Ames" = id.ames, "Goodman" = id.goodman)) }
#   }
# }
# df.low.ibs # there are four lines with IBS < 0.8 (note that all other pairs have IBS > 0.9, meaning that these four have extremely low IBS)

