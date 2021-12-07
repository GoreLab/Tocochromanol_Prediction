library(data.table)
map.a <- fread("RESULT/1.7-UseFeatures/map_a.csv")
map.b <- fread("RESULT/1.7-UseFeatures/map_b.csv")
tmp <- cbind(map.a, map.b)
df.save <- data.frame("chr" = tmp$chr,
                      "pos" = tmp$pos,
                      "MAF" = tmp$maf,
                      "GERP" = tmp$gerp_non_negative_score,
                      "MNASE_SHOOT" = tmp$mnase_hotspot_shoots,
                      "MNASE_ROOT" = tmp$mnase_hotspot_roots,
                      "REC.RATE" = tmp$r,
                      "PROX.1Kb" = tmp$gene.proximity.1Kb,
                      "PROX.10Kb" = tmp$gene.proximity.10Kb,
                      "PROX.100Kb" = tmp$gene.proximity.100Kb)
fwrite(df.save, "RESULT/1.7-UseFeatures/map_with_features.csv")
