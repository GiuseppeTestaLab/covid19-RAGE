
setwd("functional_enrichment_analysis")
source("FunctEnrich_functions.R")


#### _____----
## Upload DEGs DataMatrix----
degs <- read.table("data/DEGS.DataMatrix.csv.gz", sep = ",", header = T)

rownames(degs) <- as.vector(as.matrix(degs[, 1]))

# _ convert dataframe to numeric matrix----
degs <- data.matrix(degs[, -1])


# _ DEGs metadata----
degs.md <- as.data.frame(t(sapply(
  colnames(degs),
  function(x) {
    x1 <- unlist(strsplit(x, "__"))
    x2 <- unlist(strsplit(x1, "_"))
    le <- length(x2)
    x3 <- c(
      paste(x2[3:(le - 4)], sep = "_", collapse = "_"),
      x1[2],
      paste(x2[(le - 3):(le - 1)], sep = "_", collapse = "_")
    )
    return(x3)
  }
)))

colnames(degs.md) <- c("cells", "stat", "trend")



#### _____----
## features selection----
# _ statistics column indices----
stat.idx <- sapply(unique(degs.md$stat), function(x) list(grep(x, degs.md$stat)))


# _ matrix indices of features selected according to FDR and logFC values----
sel_feat.idx <- which(
  degs[, stat.idx$FDR] <= 5e-2 & 
    abs(degs[, stat.idx$logFC]) >= log2(2) 
)


# _ selected features indices in the logFC DataMatrix----
sel_feat.mat.idx <- as.data.frame(mat.idx(degs[, stat.idx$logFC], sel_feat.idx))

sel_feat.mat.idx <- cbind(sel_feat.mat.idx, signFC = rep("", nrow(sel_feat.mat.idx)))

sel_feat.mat.idx$signFC[which(sel_feat.mat.idx$value < 0)] <- "negFC"

sel_feat.mat.idx$signFC[which(sel_feat.mat.idx$value > 0)] <- "posFC"


# _ selected features indices in the DataMatrix----
feats.mat.idx <- match(
  colnames(degs[, stat.idx$logFC])[sel_feat.mat.idx$col],
  colnames(degs)
)


# _ features indices ----
sel_feat.mat.idx <- cbind(sel_feat.mat.idx, col.degs = feats.mat.idx)


# _ frequencies of unique features column indices ----
sel_feat.mat.idx.cu <- unique(sel_feat.mat.idx$col.degs)



####_____----
## stringdb.enrich----
# _ genelists by CELL_FAMILIES----
degs.cn.cells <- sapply(colnames(degs), function(x) strsplit(x, "__")[[1]][1])[sel_feat.mat.idx.cu]
names(degs.cn.cells) <- NULL

genelists.micropop.FC <- sapply(
  1:length(sel_feat.mat.idx.cu),
  function(x) { # x=
    idx.x <- which(sel_feat.mat.idx$col.degs %in% sel_feat.mat.idx.cu[x])
    if (length(idx.x) > 0) x1 <- rownames(degs)[sel_feat.mat.idx$row[idx.x]] else x1 <- ""
    idx.x1 <- which(!is.na(x1))

    x1a <- sapply(
      x1,
      function(y) {
        y1 <- strsplit(y, "\\|")[[1]]
        if (length(y1) > 1) y2 <- y1[2] else y2 <- y1[1]
        return(y2)
      }
    )

    x2 <- x1a[idx.x1]
    names(x2) <- NULL

    x2b <- rep(degs.cn.cells[x], length(x2))

    idx.x2 <- which(sel_feat.mat.idx$signFC[idx.x[idx.x1]] %in% "negFC")

    if (length(idx.x2) > 0) x2b[idx.x2] <- gsub("_1_2_3", "_3_2_1", x2b[idx.x2])

    x3 <- split(x2, as.factor(x2b))

    return(list(x3))
  }
)


## _ rename list elements for better readability----
genelists.micropop <- unlist(genelists.micropop.FC, recursive = F)

names(genelists.micropop) <- gsub("cell_families_", "", names(genelists.micropop))
names(genelists.micropop) <- gsub("_cells", "", names(genelists.micropop))
names(genelists.micropop) <- gsub("mild_", "Mild_", names(genelists.micropop))
names(genelists.micropop) <- gsub("s_c_", "SevCrit_", names(genelists.micropop))
names(genelists.micropop) <- gsub("_1_2_3", "__123", names(genelists.micropop))
names(genelists.micropop) <- gsub("_3_2_1", "__321", names(genelists.micropop))

	
# _ perform enrichment analysis on each cell population genelist----
enrich.li <- sapply(1:length(genelists.micropop), function(x) {
  # __ mapping gene symbols to STRING gene identifiers----
  hits <- string_db$map(
    as.data.frame(genelists.micropop[x]),
    names(genelists.micropop)[x],
    removeUnmappedRows = T
  )$STRING_id

  # __ enrichment for each of the genesets categories using STRING identifiers----
  if (length(hits) > 0) {
    x1 <- sapply(
      c("Process", "Component", "Function", "KEGG"),
      function(y) {
        y1 <- as.matrix(string_db$get_enrichment(hits, category = y))
        if (nrow(y1) > 0) y2 <- y1 else y2 <- c()
        return(y2)
      }
    )

    # __ index elements with non-null enrichment results----
    idx.x1 <- which(unlist(lapply(x1, function(z) !is.null(z))))

    if (length(idx.x1) > 0) {
      x2 <- x1[idx.x1]
      x3 <- do.call(rbind, lapply(x2, matrix, ncol = ncol(x2[[1]]), byrow = F))
      colnames(x3) <- colnames(x2[[1]])
      x3 <- as.data.frame(x3)

      # __ order merged enrichment results by FDR----
      if (nrow(x3) > 1) x3 <- x3[order(as.numeric(x3$fdr)), ]
      
      return(list(x3))
    } else {
      list(c())
    }
  } else {
    list(c())
  }
})

names(enrich.li) <- names(genelists.micropop)


# _ subset non-empty enrich.li elements----
enrich.li2 <- enrich.li[which(
  lapply(
    enrich.li,
    function(x) length(which(!is.null(dim(x))))
  ) > 0
)]


enrich.li3 <- sapply(1:length(enrich.li2), function(x) { 
  x1 <- enrich.li2[[x]]
  x2 <- cbind(rep(names(enrich.li2)[x], nrow(x1)), x1)
  colnames(x2)[1] <- "ID"
  x3 <- as.matrix(x2)

  return(x3)
})


# - enrichment results list to dataframe----
enrich.dm <- as.data.frame(
  do.call(
    rbind,
    lapply(
      enrich.li3,
      matrix,
      ncol = ncol(enrich.li3[[1]]),
      byrow = F
    )
  )
)


# - set as.numeric columns w/ numeric values----
colnames(enrich.dm) <- colnames(enrich.li3[[1]])

enrich.dm$number_of_genes <- as.numeric(enrich.dm$number_of_genes)
enrich.dm$number_of_genes_in_background <- as.numeric(enrich.dm$number_of_genes_in_background)
enrich.dm$p_value <- as.numeric(enrich.dm$p_value)
enrich.dm$fdr <- as.numeric(enrich.dm$fdr)


# - reorder enrichment results dataframe----
enrich.dm <- enrich.dm[order(
  enrich.dm$ID,
  enrich.dm$fdr,
  enrich.dm$p_value,
  -enrich.dm$number_of_genes,
  -enrich.dm$number_of_genes_in_background
), ]



####_____----
# subset entries in enrichment results by FDR----
enrich.dm2 <- enrich.dm[which(enrich.dm$fdr <= 0.001),]


# save to file selected enrichment results----
write.csv(enrich.dm2, "data/SupplementaryTable5.csv", row.names = F)



