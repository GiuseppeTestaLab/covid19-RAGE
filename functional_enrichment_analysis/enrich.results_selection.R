
setwd("functional_enrichment_analysis")
source("Covidiamo_functions.R")


# upload functional enrichment results----
funct.enrich <- read.table("data/funct.enrich.tsv.gz", sep = "\t", header = T, quote='')


# select entries in funct.enrich by FDR, max number of background genes and min number of genes per enriched functional term----
idx.enrich.fdr <- which(as.numeric(as.vector(funct.enrich$fdr)) <= 0.001)


# select entries in funct.enrich by single trend patterns----
idx.enrich.single_patterns <- setdiff(
  1:nrow(funct.enrich),
  unique(unlist(sapply(
    c("1\\.1", "1\\.2", "1\\.3", "2\\.1", "2\\.2", "2\\.3", "3\\.1", "3\\.2", "3\\.3"),
    function(x) {
      grep(x, funct.enrich$ID)
    }
  )))
)


# select entries in funct.enrich by "Mild" or "SevCrit" severity----
idx.enrich.sev <- unique(unlist(
  sapply(
    c("Mild", "SevCrit"),
    function(x) {
      grep(x, funct.enrich[, 1])
    }
  )
))


# select entries in funct.enrich by increasing or decreasing trend of expression----
idx.enrich.trend <- unique(unlist(
  sapply(
    c("123", "321"),
    function(x) {
      grep(x, funct.enrich[, 1])
    }
  )
))


# select entries in funct.enrich intersecting all the above criteria----
idx.enrich <- Reduce(
  "intersect",
  lapply(
    list(
      idx.enrich.fdr,
      idx.enrich.single_patterns,
      idx.enrich.sev,
      idx.enrich.trend
    ),
    function(x) x
  )
)


# merge GO IDs and GO description----
termID <- apply(funct.enrich[, c(3, 11)], 1, paste, sep = ": ", collapse = ": ")


# filtered funct.enrich----
funct.enrich.02 <- funct.enrich[idx.enrich, ]


# save selected enrichment table to Supplementary Table 5----
write.csv(funct.enrich.02, "data/SupplementaryTable5.csv", row.names = F)


