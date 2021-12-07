# library
library(data.table)
library(gdata)

# mkdir 
dir.save <- "RESULT/7.1-MakeExprData"
dir.create(dir.save, recursive = TRUE)

# information for the phenotype data
GbInfo <- fread("RESULT/1.2-MakeGbPhenoData/GbInfo.csv", data.table = F)
head(GbInfo)

#
f <- "RAWDATA/ExprData_282/kremling_kern_RLOG_count_matrix_all_info.txt"
ExprData <- fread(f, data.table = F)
ExprMat <- t(as.matrix(ExprData[, 6:ncol(ExprData)]))
colnames(ExprMat) <- ExprData$gene_id

# #
# vte4.01 <- as.numeric(ExprData[ExprData$gene_id == "Zm00001d017746", 6:ncol(ExprData)])
# vte4.02 <- ExprMat[, "Zm00001d017746"]
# plot(vte4.01, vte4.02)

#
f <- "RAWDATA/ExprData_282/kremling_kernel_statistics_25Mar20.xlsx"
ExprDataInfo.v1 <- read.xls(f)

# ---------------------------------------------------------------------------- #
#
f <- "RAWDATA/ExprData_282/41586_2018_BFnature25966_MOESM2_ESM.xls"
ExprDataInfo.v2 <- read.xls(f)
ExprDataInfo.v2 <- ExprDataInfo.v2[ExprDataInfo.v2$Tissue == "Kern", -ncol(ExprDataInfo.v2)] # only kernel samples

# unique names?
length(unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name)) # 231
length(unique(ExprDataInfo.v2$HMP32Name)) # 229
length(unique(ExprDataInfo.v2$HMP32Name.modif)) # 229
length(unique(ExprDataInfo.v2$AmesName)) # 230
length(unique(ExprDataInfo.v2$PhenotypeNames)) # 231

# #
# ExprDataInfo.v2$Orig_RNA_Taxa_Name[duplicated(ExprDataInfo.v2$Orig_RNA_Taxa_Name)]
# ExprDataInfo.v2$HMP32Name[duplicated(ExprDataInfo.v2$HMP32Name)]
# 
# ExprDataInfo.v2[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "B73", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "B73Htrhm", ]
# 
# ExprDataInfo.v2[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "W22", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "W22R-rstd", ]
# 
# TAB <- table(ExprDataInfo.v2$AmesName, ExprDataInfo.v2$Orig_RNA_Taxa_Name)
# v <- apply(TAB == 0, 1, sum)
# v[v == 229]
# ExprDataInfo.v2[ExprDataInfo.v2$AmesName == "Ames27101", ]

# I decided to mainly use "Orig_RNA_Taxa_Name", because this corresponds to the Lipka's name
name.ExprDataInfo.v2 <- as.character(unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name))

# attach new name column -- this will be used to match with other data sets
ExprDataInfo.v2$Name.New <- ExprDataInfo.v2$Orig_RNA_Taxa_Name
ExprDataInfo.v2$Name.New.comment <- NA

# At first, I will use the complete match
name.Lipka <- GbInfo$Name.Lipka
ExprDataInfo.v2$Name.New[!(ExprDataInfo.v2$Name.New %in% name.Lipka)] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Name.New %in% name.Lipka] <- "Perfect Match"

# manually curate the remaining...
ExprDataInfo.v2$Orig_RNA_Taxa_Name[is.na(ExprDataInfo.v2$Name.New)] # visualize
name.Lipka[!(name.Lipka %in% ExprDataInfo.v2$Name.New)] # visualize

# (1)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ab28A"] <- "AB28A"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ab28A"] <- "uppercase"

# (2)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Hi27"] <- "HI27"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Hi27"] <- "uppercase"

# (3)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky226"] <- "KY226"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky226"] <- "uppercase"

# (4)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky228"] <- "KY228"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky228"] <- "uppercase"

# (5)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki21"] <- "KI21"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki21"] <- "uppercase"

# (6)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo18W"] <- "MO18W"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo18W"] <- "uppercase"

# (7)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo45"] <- "MO45"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo45"] <- "uppercase"

# (8)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo47"] <- "MO47"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo47"] <- "uppercase"

# (9)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mp339"] <- "MP339"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mp339"] <- "uppercase"

# (10)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mt42"] <- "MT42"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mt42"] <- "uppercase"

# (11)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "N28Ht"] <- "N28HT"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "N28Ht"] <- "uppercase"

# (12)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh603"] <- "OH603"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh603"] <- "uppercase"

# (13)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi16"] <- "TZI16"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi16"] <- "uppercase"

# (14)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi25"] <- "TZI25"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi25"] <- "uppercase"

# (15,16) there are two "VaW6"
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "VaW6"] <- "VAW6"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "VaW6"] <- "uppercase"

# (17)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Yu796NS"] <- "YU796_NS"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Yu796NS"] <- "uppercase and underscore"

# (18) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR19_Santo_Domingo"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR19_Santo_Domingo"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR19_Santo_Domingo"]

# (19) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR20_Shoepeg"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR20_Shoepeg"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "MR20_Shoepeg"]

# (20) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML330"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML330"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML330"]

# (21) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML411"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML411"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML411"]

# (22) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML418"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML418"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML418"]

# (23) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML505"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML505"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML505"]

# (24) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML84"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML84"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML84"]

# (25) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML85"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML85"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML85"]

# (26) this is not included in the Lipka's phenotype file (and not found in the 282 paper)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML96"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML96"] <- ExprDataInfo.v2$Subpopulation[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML96"]

# (27) this is included in the Lipka's phenotype file but not listed in the genotype file
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CI7"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CI7"] <- "not found (no genotype)"

# (28)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky21"] <- "KY21"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ky21"] <- "uppercase"

# (29) this is not included in the Lipka's file (though included in the 282 list)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "M162W"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "M162W"] <- "not found"

# (30)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo1W"] <- "MO1W"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo1W"] <- "uppercase"

# (31)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo17"] <- "MO17"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo17"] <- "uppercase"

# (32)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo44"] <- "MO44"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo44"] <- "uppercase"

# (33)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo.G"] <- "MOG"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Mo.G"] <- "uppercase and remove dot"

# (34)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh7B"] <- "OH7B"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh7B"] <- "uppercase"

# (35)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh40B"] <- "OH40B"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh40B"] <- "uppercase"

# (36)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh43E"] <- "OH43E"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh43E"] <- "uppercase"

# (37)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa875"] <- "PA875"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa875"] <- "uppercase"

# (38)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa880"] <- "PA880"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa880"] <- "uppercase"

# (39)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa91"] <- "PA91"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Pa91"] <- "uppercase"

# (40)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va102"] <- "VA102"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va102"] <- "uppercase"

# (41)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va14"] <- "VA14"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va14"] <- "uppercase"

# (42)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va17"] <- "VA17"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va17"] <- "uppercase"

# (43)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va22"] <- "VA22"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va22"] <- "uppercase"

# (44, 45)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va35"] <- "VA35"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va35"] <- "uppercase"

# (46)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va59"] <- "VA59"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va59"] <- "uppercase"

# (47, 48)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va85"] <- "VA85"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va85"] <- "uppercase"

# (49)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va99"] <- "VA99"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Va99"] <- "uppercase"

# (50)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Wf9"] <- "WF9"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Wf9"] <- "uppercase"

# (51)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh43"] <- "OH43"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Oh43"] <- "uppercase"

# (52)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "W22R-rstd"] <- "W22_R-r:std"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "W22R-rstd"] <- "add colon"

# (53) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "4722"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "4722"] <- "popcorn"

# (54) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "HP301"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "HP301"] <- "popcorn"

# (55) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "I29"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "I29"] <- "popcorn"

# (56) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS28"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS28"] <- "popcorn"

# (57) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS69"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS69"] <- "popcorn"

# (58) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS91"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IDS91"] <- "popcorn"

# (59) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "SA24"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "SA24"] <- "popcorn"

# (60) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Sg1533"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Sg1533"] <- "popcorn"

# (61) this is included in the Lipka's phenotype file but not listed in the genotype file (popcorn)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "SG18"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "SG18"] <- "popcorn"

# (62)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "B73Htrhm"] <- "B73HTRHM"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "B73Htrhm"] <- "uppercase"

# (63) this is sweetcorn
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IA2132"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "IA2132"] <- "sweet"

# (64) maybe IL101? anyway it is sweetcorn
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il101T"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il101T"] <- "sweet"

# (65, 66) this is sweetcorn
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il14H"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il14H"] <- "sweet"

# (67) this is sweetcorn
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "P39"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "P39"] <- "sweet"

# (68) this is sweetcorn
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il677a"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Il677a"] <- "sweet"

# (69, 70) this is not included in the Lipka's file (though included in the 282 list)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML228"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML228"] <- "not found"

# (71) this is not included in the Lipka's file (though included in the 282 list)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML69"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "CML69"] <- "not found"

# (72)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki44"] <- "KI44"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki44"] <- "uppercase"

# (73) this is not included in the Lipka's file (though included in the 282 list)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "NC304"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "NC304"] <- "not found"

# (74)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi11"] <- "TZI11"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi11"] <- "uppercase"

# (75)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi18"] <- "TZI18"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Tzi18"] <- "uppercase"

# (76) this is not included in the Lipka's file (though included in the 282 list)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "A272"] <- NA
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "A272"] <- "not found"

# (77)
ExprDataInfo.v2$Name.New[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki3"] <- "KI3"
ExprDataInfo.v2$Name.New.comment[ExprDataInfo.v2$Orig_RNA_Taxa_Name == "Ki3"] <- "uppercase"

# table(table(ExprDataInfo.v2$Orig_RNA_Taxa_Name))
# table(table(ExprDataInfo.v2$Orig_RNA_Taxa_Name))
# table(table(ExprDataInfo.v2$Orig_RNA_Taxa_Name[is.na(ExprDataInfo.v2$Name.New)]))
# table(table(ExprDataInfo.v2$Orig_RNA_Taxa_Name[!is.na(ExprDataInfo.v2$Name.New)]))

# length(unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name))
# length(unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name[is.na(ExprDataInfo.v2$Name.New)]))
# length(unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name[!is.na(ExprDataInfo.v2$Name.New)]))

# ExprDataInfo.v2[ExprDataInfo.v2$Name.New.comment == "landrace_or_non282_CIMMYT_line", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Name.New.comment == "popcorn", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Name.New.comment == "sweet", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Name.New.comment == "not found", ]
# ExprDataInfo.v2[ExprDataInfo.v2$Name.New.comment == "not found (no genotype)", ]

# ---------------------------------------------------------------------------- #
# make 282 expression data
a <- as.character(unique(ExprDataInfo.v2$Name.New[!is.na(ExprDataInfo.v2$Name.New)]))
GbExprInfo <- NULL
for ( i in 1:length(a) ) {
	name.i <- a[i]
	tf.i <- ExprDataInfo.v2$Name.New %in% name.i
	df.i <- data.frame("Name.Lipka.Manual" = name.i, 
										 "Orig.RNA.Taxa.Name.Nature.2018" = ExprDataInfo.v2$Orig_RNA_Taxa_Name[tf.i],
										 "Orig.RNA.Sample.Name.Nature.2018" = ExprDataInfo.v2$OrigColNameFromRNAExpressionValue[tf.i],
										 "Reads.per.file.Nature.2018" = ExprDataInfo.v2$Reads.per.file[tf.i])
	df.i$Sample.Name.Josh <- paste0(df.i$Orig.RNA.Taxa.Name.Nature.2018, "_Kern")
	if ( "Mo.G_Kern" %in% df.i$Sample.Name.Josh ) { df.i$Sample.Name.Josh <- "Mo_G_Kern" } # exception handling
	ExprDataInfo.v1.i <- ExprDataInfo.v1[ExprDataInfo.v1$Sample.Name. %in% unique(df.i$Sample.Name.Josh), ]
	n.reads.josh <- as.integer(gsub(",", "", (ExprDataInfo.v1.i$Initial.Number.of.Reads)))
	m <- match(df.i$Reads.per.file.Nature.2018, n.reads.josh)
	RNA.initial.n.read.Josh <- as.integer(gsub(",", "", (ExprDataInfo.v1.i$Initial.Number.of.Reads[m])))
	RNA.clean.n.read.Josh <- as.integer(gsub(",", "", (ExprDataInfo.v1.i$Number.Cleaned.Reads[m])))
	RNA.Sample.Name.Josh <- ExprDataInfo.v1.i$Mapping.Name..SRA.ID...Pedigree.[m]
	Overall.Alignment.Rate.Josh <- as.numeric(gsub("%", "", ExprDataInfo.v1.i$Overall.Alignment.Rate.[m])) / 100
	df.i$Sample.Name.Full.Josh <- RNA.Sample.Name.Josh
	df.i$Initial.Number.of.Reads.Josh <- RNA.initial.n.read.Josh
	df.i$Number.Cleaned.Reads.Josh <- RNA.clean.n.read.Josh
	df.i$Overall.Alignment.Rate.Josh <- Overall.Alignment.Rate.Josh
	df.i$Comment.Match.Name.from.Nature.to.Lipka <- ExprDataInfo.v2$Name.New.comment[tf.i]
	GbExprInfo <- rbind.data.frame(GbExprInfo, df.i)
}
GbExprInfo.fin <- GbExprInfo[, c(3, 2, 4, 5:9, 1, 10)]

# check uniqueness
length(GbExprInfo.fin$Sample.Name.Full.Josh) # 223
length(unique(GbExprInfo.fin$Sample.Name.Full.Josh)) # 223

# save
f <- paste0(dir.save, "/GbExprDataInfo_full.csv")
write.csv(GbExprInfo.fin, f, row.names = F)

# # Summary
# table(GbExprInfo.fin$Orig.RNA.Taxa.Name.Nature.2018)
# table(GbExprInfo.fin$Comment.Match.Name.from.Nature.to.Lipka)
# length(unique(GbExprInfo.fin$Orig.RNA.Taxa.Name.Nature.2018[GbExprInfo.fin$Comment.Match.Name.from.Nature.to.Lipka == "Perfect Match"]))
# length(unique(GbExprInfo.fin$Orig.RNA.Taxa.Name.Nature.2018[GbExprInfo.fin$Comment.Match.Name.from.Nature.to.Lipka != "Perfect Match"]))


# ---------------------------------------------------------------------------- #
# Retain the best quality one
name.lipka.all <- as.character(unique(GbExprInfo.fin$Name.Lipka.Manual))
GbExprInfo.unique <- NULL
for ( i in 1:length(name.lipka.all) ) {
	name.i <- name.lipka.all[i]
	tf.i <- GbExprInfo.fin$Name.Lipka.Manual == name.i
	max.num.i <- which.max(GbExprInfo.fin$Overall.Alignment.Rate.Josh[tf.i])
	df.i <- GbExprInfo.fin[tf.i, ][max.num.i, ]
	GbExprInfo.unique <- rbind.data.frame(GbExprInfo.unique, df.i)
}
nrow(GbExprInfo.unique) # 202 unique genotype

# write
f <- paste0(dir.save, "/GbExprDataInfo_unique.csv")
write.csv(GbExprInfo.unique, f, row.names = F)



# ---------------------------------------------------------------------------- #
# correlation among duplicated samples
name.dup <- unique(ExprDataInfo.v2$Orig_RNA_Taxa_Name[duplicated(ExprDataInfo.v2$Orig_RNA_Taxa_Name)])
name.dup.2 <- paste0(name.dup, "_Kern")
avg.r.all <- c()
for ( i in 1:length(name.dup.2) ) {
	num <- grep(name.dup.2[i], rownames(ExprMat))
	cor.mat <- cor(t(ExprMat[num, ]))
	r.vec <- cor.mat[upper.tri(cor.mat)]
	avg.r.all <- c(avg.r.all, mean(r.vec))
}
df.avg.cor <- data.frame("Accession.Name" = name.dup,
												 "Correlation" = round(avg.r.all, 3))
# df.avg.cor


# ---------------------------------------------------------------------------- #
# expression data
GbExprInfo.unique$Sample.Name.Full.Josh.v2 <- paste0(GbExprInfo.unique$Sample.Name.Full.Josh, "_Kern")
m <- match(GbExprInfo.unique$Sample.Name.Full.Josh.v2, rownames(ExprMat))
ExprMat.unique <- ExprMat[m, ] # sort
all(rownames(ExprMat.unique) == GbExprInfo.unique$Sample.Name.Full.Josh.v2) # OK
ExprDat.unique <- cbind("SampleName.Josh" = rownames(ExprMat.unique),
												"Name.Lipka" = GbExprInfo.unique$Name.Lipka.Manual,
												as.data.frame(ExprMat.unique))
rownames(ExprDat.unique) <- NULL

# save
f <- paste0(dir.save, "/GbExprData.csv")
fwrite(ExprDat.unique, f)

