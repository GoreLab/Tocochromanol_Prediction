

# Function: 250Kb or 1Mb from peak SNP
MyFun.GetSnp.QTL <- function(trait, bp.width.qtl, GenData = GenData, map = map, ...) {
   # object to write
   SnpNumList <- list()
   
   # ---------------------------- 1. Make datasets ---------------------------- #
   # knwon QTLs
   GenData.trait <- GenData[GenData$Trait %in% trait, ]
   GenData.trait <- GenData.trait[!is.na(GenData.trait$Peak.Pos.v4), ] # remove SNPs which I could not get v4 position
   
   # ---------------------------- 2. Get SNPs ---------------------------- #
   # Get SNPs
   SnpNumList <- list()
   for ( j in 1:nrow(GenData.trait) ) {
      GenData.trait.j <- GenData.trait[j, ]
      name.j <- GenData.trait.j$NAME
      chr <- as.numeric(gsub("Chr" , "",GenData.trait.j$Chr))
      pos <- GenData.trait.j$Peak.Pos.v4
      bp.start <- pos - bp.width.qtl
      bp.end <- pos + bp.width.qtl
      tf <- (map$chr == chr) & (bp.start < map$pos) & (map$pos < bp.end)
      SnpNumList[[name.j]] <- which(tf)
   } # here are the Qtls
   SnpNumList[["NonQtl"]] <- setdiff(1:nrow(map), unique(unlist(SnpNumList))) # and all Qtls
   
   # return
   return(SnpNumList)
}

# Function: from genes
MyFun.GetSnp.GENE <- function(trait, bp.width.gene, TocoGene = TocoGene, map = map, ...) {
   # object to write
   SnpNumList <- list()
   
   # ---------------------------- 1. Make datasets ---------------------------- #
   # gene data of the trait
   TocoGene.tmp <- TocoGene[!is.na(TocoGene[[trait]]), ]
   TocoGene.trait <- TocoGene.tmp[ c("qtl.id", "Gene", "Chr", "Start", "End")]
   
   # ---------------------------- 2. Get SNPs ---------------------------- #
   # Get SNPs: genes
   for ( k in 1:nrow(TocoGene.trait) ) {
      TocoGene.trait.k <- TocoGene.trait[k, ]
      name.k <- TocoGene.trait.k$Gene
      chr <- TocoGene.trait.k$Chr
      bp.start <- TocoGene.trait.k$Start - bp.width.gene
      bp.end <- TocoGene.trait.k$End + bp.width.gene
      tf <- (map$chr == chr) & (bp.start < map$pos) & (map$pos < bp.end)
      if ( sum(tf) != 0 ) {
         SnpNumList[[name.k]] <- which(tf)
      } 
   }
   
   # remainings
   SnpNumList[["NonQtl"]] <- setdiff(1:nrow(map), unique(unlist(SnpNumList))) # and all Qtls
   
   # return
   return(SnpNumList)
}

# Function: support interval
MyFun.GetSnp.SI <- function(trait, GenData = GenData, map = map, ...) {
   # object to write
   SnpNumList <- list()
   
   # ---------------------------- 1. Make datasets ---------------------------- #
   # knwon QTLs
   GenData.trait <- GenData[GenData$Trait %in% trait, ]
   tf.L <- is.na(GenData.trait$SupInt.L.Pos.v4)
   tf.R <- is.na(GenData.trait$SupInt.R.Pos.v4)
   tf <- !(tf.L | tf.R)
   GenData.trait <- GenData.trait[tf, ] # remove SNPs which I could not get SI
   
   # ---------------------------- 2. Get SNPs ---------------------------- #
   # Get SNPs
   SnpNumList <- list()
   for ( j in 1:nrow(GenData.trait) ) {
      GenData.trait.j <- GenData.trait[j, ]
      name.j <- GenData.trait.j$NAME
      chr <- as.numeric(gsub("Chr" , "",GenData.trait.j$Chr))
      bp.start <- GenData.trait.j$SupInt.L.Pos.v4
      bp.end <- GenData.trait.j$SupInt.R.Pos.v4
      tf <- (map$chr == chr) & (bp.start < map$pos) & (map$pos < bp.end)
      if ( sum(tf) != 0 ) {
         SnpNumList[[name.j]] <- which(tf)
      } 
   } # here are the Qtls
   SnpNumList[["NonQtl"]] <- setdiff(1:nrow(map), unique(unlist(SnpNumList))) # and all Qtls
   
   # return
   return(SnpNumList)
}

# Function: hybrid
MyFun.GetSnp.HYBRID <- function(trait, bp.width.qtl, bp.width.gene,
                                TocoGene = TocoGene, GenData = GenData, map = map, ...) {
   # object to write
   SnpNumList <- list()
   
   # ---------------------------- 1. Make datasets ---------------------------- #
   # gene data of the trait
   TocoGene.tmp <- TocoGene[!is.na(TocoGene[[trait]]), ]
   TocoGene.trait <- TocoGene.tmp[ c("qtl.id", "Gene", "Chr", "Start", "End")]
   
   # QTL data of the trait
   GenData.trait <- GenData[GenData$Trait %in% trait, ]
   GenData.trait <- GenData.trait[!is.na(GenData.trait$Peak.Pos.v4), ] # remove SNPs which I could not get v4 position
   
   # remove QTL of the known genes
   GenData.trait$QTL.ID <- as.character(GenData.trait$QTL.ID)
   id.qtl <- setdiff(GenData.trait$QTL.ID, TocoGene.trait$qtl.id)
   GenData.trait <- GenData.trait[GenData.trait$QTL.ID %in% id.qtl, ]
   
   # ---------------------------- 2. Get SNPs ---------------------------- #
   # Get SNPs: QTLs
   SnpNumList <- list()
   for ( j in 1:nrow(GenData.trait) ) {
      GenData.trait.j <- GenData.trait[j, ]
      name.j <- GenData.trait.j$NAME
      chr <- as.numeric(gsub("Chr" , "",GenData.trait.j$Chr))
      pos <- GenData.trait.j$Peak.Pos.v4
      bp.start <- pos - bp.width.qtl
      bp.end <- pos + bp.width.qtl
      tf <- (map$chr == chr) & (bp.start < map$pos) & (map$pos < bp.end)
      if ( sum(tf) != 0 ) {
         SnpNumList[[name.j]] <- which(tf)
      }
   }
   
   # Get SNPs: genes
   for ( k in 1:nrow(TocoGene.trait) ) {
      TocoGene.trait.k <- TocoGene.trait[k, ]
      name.k <- TocoGene.trait.k$Gene
      chr <- TocoGene.trait.k$Chr
      bp.start <- TocoGene.trait.k$Start - bp.width.gene 
      bp.end <- TocoGene.trait.k$End + bp.width.gene
      tf <- (map$chr == chr) & (bp.start < map$pos) & (map$pos < bp.end)
      if ( sum(tf) != 0 ) {
         SnpNumList[[name.k]] <- which(tf)
      } 
   }
   
   # remainings
   SnpNumList[["NonQtl"]] <- setdiff(1:nrow(map), unique(unlist(SnpNumList))) # and all Qtls
   
   # return
   return(SnpNumList)
}






