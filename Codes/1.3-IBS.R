library(openxlsx)

# data frame of 1462 Ames lines (from Di's GWAS/TWAS paper)
df.GwasTwasSuppl <- openxlsx::read.xlsx("Supplementary Table S1.xlsx", startRow = 2)
df.GwasTwasSuppl <- df.GwasTwasSuppl[df.GwasTwasSuppl$Included.in.GWAS == "YES", c("Accession.number", "Inbred.line.name")]

# data frame of my phenotype data 
df.Good <- read.csv("GbPheno.csv")
df.Ames <- read.csv("AmesPheno.csv")

# check
setequal(df.Ames$ID, df.GwasTwasSuppl$Accession.number) # OK
length(df.Ames$ID) == length(unique(df.Ames$ID)) # OK

# Ames sample info
df.Ames <- data.frame("Accession.Number" = df.Ames$ID,
                      "Inbred.Line.Name.Ames" = df.GwasTwasSuppl$Inbred.line.name[match(df.Ames$ID, df.GwasTwasSuppl$Accession.number)])

# Goodman sample info
df.Good <- data.frame("GBS.Sample" = df.Good$ID,
                      "Inbred.Line.Name.Goodman" = substr(df.Good$ID, 1, nchar(df.Good$ID) - 10))
df.Good$Inbred.Line.Name.Goodman[df.Good$Inbred.Line.Name.Goodman == "MS71"] <- "Ms71" # manual correction

# write the two tables in csv file -- and manually go through all Goodman names
dir.create("SampleInfo")
write.table(df.Ames[order(df.Ames$Inbred.Line.Name.Ames), ], "SampleInfo/AmesSamples.txt", row.names = F, quote = F)
write.table(df.Good, "SampleInfo/GoodmanSamples.txt", row.names = F, quote = F)

# ---------------------------------------------------------------------------- #
# ----- Do manual comparison using the files from above
# ---------------------------------------------------------------------------- #
# ...
# ...
# ---------------------------------------------------------------------------- #
# ----- manual comparison was done.
# ---------------------------------------------------------------------------- #


dir.create("OverlapResult")
# load the result of my manual comparison
df.match <- read.xlsx("SampleInfo/MatchSamples.xlsx")

# 
library(data.table)
df.geno.012 <- fread("merged_ames_282_LD0.1_all_chr.012", data.table = F)
df.indv <- fread("merged_ames_282_LD0.1_all_chr.012.indv", data.table = F, header = F)

# df to mat
mat.geno.012 <- as.matrix(df.geno.012[, -1])
rownames(mat.geno.012) <- df.indv$V1
P <- ncol(mat.geno.012)

# 
df.match$IBS.02 <- df.match$IBS.01 <- NA
for ( i in 1:nrow(df.match) ) {
  flag.i <- df.match$Overlap.Manually.Found[i]
  good.gbs.i <- df.match$GBS.Sample[i]
  num.gbs.i <- which(df.indv$V1 == good.gbs.i)
  
  if ( flag.i == "Yes" ) {
    # calc IBS
    ames.name.01.i <- df.match$Manually.Matched.Ames.Line.Name[i]
    ames.gbs.01.i <- df.match$Manually.Identified.Ames.Accession.Name[i]
    checker <- df.Ames$Inbred.Line.Name.Ames[df.Ames$Accession.Number == ames.gbs.01.i]
    if (checker != ames.name.01.i) { print(i) }
    num.ames.01.i <- which(df.indv$V1 == ames.gbs.01.i)
    snp.good.i <- mat.geno.012[num.gbs.i, ]
    snp.ames.01.i <- mat.geno.012[num.ames.01.i, ]
    ibs.01.i <- sum(snp.good.i == snp.ames.01.i) / P
    df.match$IBS.01[i] <- ibs.01.i
    
    # calc IBS for the 2nd candidate if there is
    ames.name.02.i <- df.match$Manually.Matched.Ames.Line.Name.02[i]
    ames.gbs.02.i <- df.match$Manually.Identified.Ames.Accession.Name.02[i]
    if ( !is.na(ames.name.02.i) ) {
      checker <- df.Ames$Inbred.Line.Name.Ames[df.Ames$Accession.Number == ames.gbs.02.i]
      if (checker != ames.name.02.i) { print(i) }
      num.ames.02.i <- which(df.indv$V1 == ames.gbs.02.i)
      snp.good.i <- mat.geno.012[num.gbs.i, ]
      snp.ames.02.i <- mat.geno.012[num.ames.02.i, ]
      ibs.02.i <- sum(snp.good.i == snp.ames.02.i) / P
      df.match$IBS.02[i] <- ibs.02.i
    }
  }
}

# Ames lines
df.Ames$Inbred.Line.Name.Goodman <- NA
df.Ames$Pairwise.IBS <- NA
df.Ames$Overlap <- NA
for ( i in 1:nrow(df.Ames) ) {
  gbs.ames.i <- df.Ames$Accession.Number[i]
  name.i <- df.Ames$Inbred.Line.Name.Ames[i]
  
  if ( !(gbs.ames.i %in% df.match$Manually.Identified.Ames.Accession.Name) & !(gbs.ames.i %in% df.match$Manually.Identified.Ames.Accession.Name.02 ) ) {
    df.Ames$Overlap[i] <- "No"
  }
  
  if ( gbs.ames.i %in% df.match$Manually.Identified.Ames.Accession.Name ) {
    m <- match(gbs.ames.i, df.match$Manually.Identified.Ames.Accession.Name)
    df.Ames$Inbred.Line.Name.Goodman[i] <- df.match[m, "Inbred.Line.Name.Goodman"]
    df.Ames$Pairwise.IBS[i] <- round(df.match[m, "IBS.01"], 3)
    if ( df.match[m, "IBS.01"] >= 0.8 ) {
      df.Ames$Overlap[i] <- "Yes"
    } else {
      df.Ames$Overlap[i] <- "No (IBS < 0.8)"
    }
  }
  
  if ( gbs.ames.i %in% df.match$Manually.Identified.Ames.Accession.Name.02 ) {
    m <- match(gbs.ames.i, df.match$Manually.Identified.Ames.Accession.Name.02)
    df.Ames$Inbred.Line.Name.Goodman[i] <- df.match[m, "Inbred.Line.Name.Goodman"]
    df.Ames$Pairwise.IBS[i] <- round(df.match[m, "IBS.02"], 3)
    if ( df.match[m, "IBS.02"] >= 0.8 ) {
      df.Ames$Overlap[i] <- "Yes"
    } else {
      df.Ames$Overlap[i] <- "No (IBS < 0.8)"
    }
  }
}
write.csv(df.Ames, "OverlapResult/AmesOverlapResult.csv", row.names = F)

# Goodman lines
df.Good$Inbred.Line.Name.Ames <- NA
df.Good$Pairwise.IBS <- NA
df.Good$Overlap <- NA
for ( i in 1:nrow(df.Good) ) {
  gbs.good.i <- df.Good$GBS.Sample[i]
  df.match.sub <- df.match[df.match$GBS.Sample == gbs.good.i, ]
  
  if ( df.match.sub$Overlap.Manually.Found == "No" ) { # no overlap
    df.Good$Overlap[i] <- "No"
  } else if ( is.na(df.match.sub$Manually.Matched.Ames.Line.Name.02) ) { # one overlap(?)
    df.Good$Inbred.Line.Name.Ames[i] <- df.match.sub$Manually.Matched.Ames.Line.Name
    df.Good$Pairwise.IBS[i] <- round(df.match.sub$IBS.01, 3)
    if ( df.match.sub$IBS.01 >= 0.8 ) {
      df.Good$Overlap[i] <- "Yes"
    } else {
      df.Good$Overlap[i] <- "No (IBS < 0.8)"
    }
  } else { # two overlaps(?)
    df.Good$Inbred.Line.Name.Ames[i] <- paste0(df.match.sub$Manually.Matched.Ames.Line.Name, "(",
                                               df.match.sub$Manually.Identified.Ames.Accession.Name, ") and ",
                                               df.match.sub$Manually.Matched.Ames.Line.Name.02, "(",
                                               df.match.sub$Manually.Identified.Ames.Accession.Name.02, ")")
    df.Good$Pairwise.IBS[i] <- paste0(round(df.match.sub$IBS.01, 3), "(",
                                      df.match.sub$Manually.Identified.Ames.Accession.Name, ") and ",
                                      round(df.match.sub$IBS.02, 3), "(",
                                      df.match.sub$Manually.Identified.Ames.Accession.Name.02, ")")
    if ( (df.match.sub$IBS.01 >= 0.8) & (df.match.sub$IBS.02 >= 0.8) ) {
      df.Good$Overlap[i] <- "Yes (with two Ames GBS samples)"
    } else if ( (df.match.sub$IBS.01 >= 0.8) & (df.match.sub$IBS.02 < 0.8) ) {
      df.Good$Overlap[i] <- paste0("Yes (with one Ames GBS sample ", df.match.sub$Manually.Identified.Ames.Accession.Name, ")")
    } else if ( (df.match.sub$IBS.01 < 0.8) & (df.match.sub$IBS.02 >= 0.8) ) {
      df.Good$Overlap[i] <- paste0("Yes (with one Ames GBS sample ", df.match.sub$Manually.Identified.Ames.Accession.Name.02, ")")
    } else {
      df.Good$Overlap[i] <- "No"
    }
  }
}
write.csv(df.Good, "OverlapResult/GoodmanOverlapResult.csv", row.names = F)

# make file for 1704
df.Good$Inbred.Line.Name.Consistent <- df.Good$Inbred.Line.Name.Goodman
df.Good$Inbred.Line.Name.Consistent[!is.na(df.Good$Inbred.Line.Name.Ames)] <- df.Good$Inbred.Line.Name.Ames[!is.na(df.Good$Inbred.Line.Name.Ames)]
df.1704 <- data.frame("Panel" = rep(c("Ames", "Goodman"), times = c(nrow(df.Ames), nrow(df.Good))),
                      "Inbred.Line.Name" = c(df.Ames$Inbred.Line.Name.Ames, df.Good$Inbred.Line.Name.Consistent),
                      "Accession.Number.or.GBS.sample" = c(df.Ames$Accession.Number, df.Good$GBS.Sample),
                      "Pairwise.IBS" = c(df.Ames$Pairwise.IBS, df.Good$Pairwise.IBS),
                      "Overlap" = c(df.Ames$Overlap, df.Good$Overlap))
write.csv(df.1704, "OverlapResult/OverlapResult_all.csv", row.names = F)

# ################################################################################
# # load the updated table S1
# df.match.with.Mike <- read.xlsx("Table S1.xlsx", startRow = 2)
# 
# # ---------------------- #
# # (1) on the Ames panel  
# # ---------------------- #
# df.match.with.Mike.Ames <- df.match.with.Mike[df.match.with.Mike$Panel == "Ames", ]
# df.match.by.Ryokei.Ames <- df.1704[df.1704$Panel == "Ames", ]
# 
# # check the order of lines  -- ok
# df.match.with.Mike.Ames[df.match.with.Mike.Ames$Inbred.line.name != df.match.by.Ryokei.Ames$Inbred.Line.Name, ]
# df.match.by.Ryokei.Ames[df.match.with.Mike.Ames$Inbred.line.name != df.match.by.Ryokei.Ames$Inbred.Line.Name, ]
# 
# # compare -- I have two newly defined overlapping lines
# table("Overlap.with.Mike" = df.match.with.Mike.Ames$Overlap.a, "Overlap.new.Ryokei" = df.match.by.Ryokei.Ames$Overlap)
# 
# # show the one case -- W22R
# df.match.by.Ryokei.Ames[(df.match.with.Mike.Ames$Overlap.a == "No") & (df.match.by.Ryokei.Ames$Overlap == "Yes"), ]
# 
# 
# # ---------------------- #
# # (2) on the Goodman panel  
# # ---------------------- #
# df.match.with.Mike.Good <- df.match.with.Mike[df.match.with.Mike$Panel == "Goodman", ]
# df.match.by.Ryokei.Good <- df.1704[df.1704$Panel == "Goodman", ]
# 
# # compare -- I have two newly defined overlapping lines
# table("Overlap.with.Mike" = df.match.with.Mike.Good$Overlap.a, "Overlap.new.Ryokei" = df.match.by.Ryokei.Good$Overlap)
# df.match.by.Ryokei.Good[df.match.with.Mike.Good$Overlap.a != df.match.by.Ryokei.Good$Overlap, ]
# 
