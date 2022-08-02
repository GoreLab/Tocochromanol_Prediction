# source
library(data.table)
library(readxl)

# mkdir 
dir.create("RESULT/1.3-CalcIbs", recursive = T)

# ---------------------------------------------------------------------------- #
# ----- Check IBS
# ---------------------------------------------------------------------------- #
#
df.ames.info <- data.frame(read_xlsx("vitamaize_accession_removal_20191113.xlsx"))

#
IBS.mat <- read.csv("RESULT/1.3-CalcIbs/IbsMat.csv")
IBS.mat <- as.matrix(IBS.mat[, -1])
rownames(IBS.mat) <- colnames(IBS.mat)
colnames(IBS.mat)[1463:1704] <- rownames(IBS.mat)[1463:1704] <- substr(rownames(IBS.mat)[1463:1704], 
                                                                       1, 
                                                                       nchar(rownames(IBS.mat)[1463:1704]) - 10)

# #
# df.LineName <- data.table::fread("Table S1.txt", data.table = F)
# colnames(df.LineName) <- c("Panel", "Inbred.line.name", "overlap")

# load table S1 (with Ames Accession ID for the ones duplicated within Ames)
df.LineName <- as.data.frame(read_excel("Table S1.xlsx"))
head(df.LineName)

#
name.ames <- df.LineName$Inbred.line.name[df.LineName$Panel == "Ames"]
name.good <- df.LineName$Inbred.line.name[df.LineName$Panel == "Goodman"]

#
df.name.ames <- data.frame("Name" = name.ames)
df.name.ames$Name.v2 <- df.name.ames$Name
df.name.ames$Name.v2 <- gsub(" ", "", df.name.ames$Name.v2)

#
df.name.good <- data.frame("Name" = name.good)
df.name.good$Name.v2 <- df.name.good$Name
df.name.good$Name.v2 <- gsub(" ", "", df.name.good$Name.v2)
df.name.good$Name.v2[df.name.good$Name.v2 == "Hy"] <- "Ill.Hy"
df.name.good$Name.v2[df.name.good$Name.v2 == "MoG"] <- "Mo.G"
df.name.good$Name.v2[df.name.good$Name.v2 == "W117Ht"] <- "W117HT"

#
overlapping.names <- intersect(df.name.ames$Name.v2, df.name.good$Name.v2)
length(overlapping.names) # 170 candidate overlaps

# Ames ID (to match GBS samples in Ames)
df.potential.overlap <- data.frame("Overlap.Name" = overlapping.names,
                                   "Name.Goodman" = df.name.good$Name[match(overlapping.names, df.name.good$Name.v2)],
                                   "Name.Ames" = df.name.ames$Name[match(overlapping.names, df.name.ames$Name.v2)],
                                   "Ames.ID.01" = NA, "Ames.ID.02" = NA)
for ( i in 1:nrow(df.potential.overlap) ) {
  name.i <- df.potential.overlap$Name.Ames[i]
  df.ames.info.i <- df.ames.info[df.ames.info$Pedigree_vita %in% name.i, ]
  if ( nrow(df.ames.info.i) == 1 ) {
    df.potential.overlap$Ames.ID.01[i] <- df.ames.info.i$Accession.Number_ion
  } else if ( nrow(df.ames.info.i) == 2 ) {
    df.potential.overlap$Ames.ID.01[i] <- df.ames.info.i$Accession.Number_ion[1]
    df.potential.overlap$Ames.ID.02[i] <- df.ames.info.i$Accession.Number_ion[2]
  }
}

# Goodman names
df.potential.overlap$Name.Goodman <- gsub(" ", "", df.potential.overlap$Name.Goodman)
df.potential.overlap$Name.Goodman <- gsub("-", ".", df.potential.overlap$Name.Goodman)
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "38.11"] <- "X38.11"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "33.16"] <- "X33.16"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "4226"] <- "X4226"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "Oh7B"] <- "OH7B"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "Ms71"] <- "MS71"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "Va102"] <- "VA102"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "DE2"] <- "DE_2"
df.potential.overlap$Name.Goodman[df.potential.overlap$Name.Goodman == "DE3"] <- "DE_3"


#
df.potential.overlap$IBS <- df.potential.overlap$IBS.02 <- df.potential.overlap$IBS.01 <- NA
for ( i in 1:nrow(df.potential.overlap) ) {
  Good.name.i <- df.potential.overlap$Name.Goodman[i]
  Ames.name.01.i <- df.potential.overlap$Ames.ID.01[i]
  Ames.name.02.i <- df.potential.overlap$Ames.ID.02[i]
  if ( is.na(Ames.name.02.i) ) {
    df.potential.overlap$IBS.01[i] <- IBS.mat[Good.name.i, Ames.name.01.i]
    df.potential.overlap$IBS[i] <- IBS.mat[Good.name.i, Ames.name.01.i]
  } else {
    df.potential.overlap$IBS.01[i] <- IBS.mat[Good.name.i, Ames.name.01.i]
    df.potential.overlap$IBS.02[i] <- IBS.mat[Good.name.i, Ames.name.02.i]
    df.potential.overlap$IBS[i] <- max(c(IBS.mat[Good.name.i, Ames.name.01.i], IBS.mat[Good.name.i, Ames.name.02.i]))
  }
}

#
hist(df.potential.overlap$IBS)
df.potential.overlap[df.potential.overlap$IBS < 0.8, ]
write.csv(df.potential.overlap, "Table_IBS_overlap_confirm.csv")


