
setwd("functional_enrichment_analysis")
source("FunctEnrich_functions.R")


#### _____----
## Upload SupplementaryTable5 DataMatrix----
enrich.dm2 <- read.csv("data/SupplementaryTable5.csv", sep = ",", header = T)


#### _____----
# further subset entries in enrichment results by min nr. of genes, max nr. of background genes and including only Gene Ontology terms----
enrich.dm3 <- enrich.dm2[
  which(
    as.numeric(as.vector(enrich.dm2$number_of_genes)) >= 3 &
      as.numeric(as.vector(enrich.dm2$number_of_genes_in_background)) <= 500 &
      !grepl("hsa", enrich.dm2$term)
  ),
]


# re-format entries (for better readability) in some fields of filtered enrich.dm3----
enrich.dm3$number_of_genes <- gsub(" ", "", enrich.dm3$number_of_genes)

enrich.dm3$number_of_genes_in_background <- gsub(" ", "", enrich.dm3$number_of_genes_in_background)


## ______----
# cell families and patterns----
enrich.dm3.famnam <-
  as.data.frame(t(
    sapply(
      enrich.dm3$ID,
      function(x) {
        x1 <- strsplit(x, "__")[[1]]
        return(x1)
      }
    )
  ))

colnames(enrich.dm3.famnam) <- c("CellFam", "Pattern")


# merge GO IDs and GO description----
termID.dm3 <- apply(enrich.dm3[, c(3, 11)], 1, paste, sep = ": ", collapse = ": ")


# frequencies table of GO terms and cells by patterns----
enrich.dm3.ta <-
  table(
    enrich.dm3$ID,
    termID.dm3
  )


# set up a data matrix w/ logged FDR values----
enrich.dm3.ta2 <- enrich.dm3.ta

enrich.dm3.ta2[which(enrich.dm3.ta2 >= 0)] <- 1


# populating the data matrix w/ logged FDR values----
for (i in rownames(enrich.dm3.ta)) {
  idx1 <- which(enrich.dm3$ID %in% i)
  idx2r <- which(rownames(enrich.dm3.ta) %in% i)
  i1 <- names(which(enrich.dm3.ta[idx2r, ] > 0))
  idx2c <- which(colnames(enrich.dm3.ta) %in% i1)
  idx3 <- which(termID.dm3[idx1] %in% i1)

  for (i2 in idx2c) {
    idx4 <- which(termID.dm3[idx1[idx3]] %in% colnames(enrich.dm3.ta)[i2])

    enrich.dm3.ta2[idx2r, i2] <- log(-log(min(
      as.numeric(enrich.dm3$fdr[idx1[idx3][idx4]]),
      na.rm = T
    ), 10), 2)
  }
}

enrich.dm3.ta3 <- as.matrix(t(enrich.dm3.ta2))


# _ reformat heatmap colnames----
heatmap.cn <- paste(
  c("Mild_", "SevCrit_"),
  gsub(
    "Myel.",
    "Myeloid",
    gsub(
      "_",
      "__",
      c("B_321", "B_321", "CP_321", "CP_321", "Myel._123", "Myel._123", "Myel._321", "Myel._321", "NK_321", "NK_321", "T_123", "T_123")
    )
  ),
  sep = ""
)


# _ select enrich.dm3.ta3 columns according to the heatmap colnames----
enrich.dm3.ta3.cn <- match(heatmap.cn,colnames(enrich.dm3.ta3))

enrich.dm3.ta4 <- enrich.dm3.ta3[,enrich.dm3.ta3.cn]

enrich.dm3.ta4[,4] <- rep(1,nrow(enrich.dm3.ta4))

colnames(enrich.dm3.ta4)[4] <- gsub("Mild","SevCrit",colnames(enrich.dm3.ta4)[3])

colnames(enrich.dm3.ta4) <- gsub("Myeloid","Myel.",colnames(enrich.dm3.ta4))

	
## ______----
# save to file selected enrichment results----
write.table(enrich.dm3.ta4, "data/Enrich_DataMatrix.tsv", sep = "\t", col.names = NA)


