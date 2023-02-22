
setwd("functional_enrichment_analysis")
source("Covidiamo_functions.R")


# upload the selected functional enrichment datamatrix----
enrich.hg <- read.table("data/enrich_datamat.tsv",sep="\t",header=T)
	
colnames(enrich.hg) <- gsub("\\.1$","",colnames(enrich.hg))


# _ define pathways families to be used in GO terms selection----
PathFam.terms <- list(
  CellCycle = c("cell cycle", "cell differe", "division", "replic", "mitotic"),
  ImmuneResponse = c("immun", "interfer", "cytokine", "nflamm", "HLA", "MHC", "virus", "leuk", "viral", "virus", "paras", "graft", "Graft", "ntigen", "llograft", "nfection", "bacter", "RAGE"),
  Endoplasm = c("endoplasmic ret", "cytoplasm", "organell", "localiz"),
  Hematopoiesis = c("ematopoie"),
  Metabolism = c("etabol", "synthe", "catabol", "anabol"),
  Mitochondrial = c("Mitoch", "mitoch", "oxphos", "xidat", "respir", "cytochrom", "proton trans", "xidoreduct"),
  Ribosome = c("Ribo", "ribo", "transl", "rRNA")
)


## _ index pathway descriptions by pathways families----
idx.PathFam <- c()
for (i in 1:length(PathFam.terms)) { 
  i1 <- as.numeric(unlist(
    sapply( 
      PathFam.terms[[i]],
      function(y) grep(y, rownames(enrich.hg))
    )
  ))
  idx.PathFam <- c(idx.PathFam, list(i1))
}

names(idx.PathFam) <- names(PathFam.terms)

idx.PathFam$Other <- setdiff(1:nrow(enrich.hg), unlist(idx.PathFam))

legend.row2 <- names(idx.PathFam)



## _ reorder the indexed pathway families descriptions for the stacked histogram----
counter <- 0
idx.PathFam.StackedHist <- c()

for (i in 1:length(idx.PathFam)) {
  i1 <- 1:length(idx.PathFam[[i]])
  i2 <- i1 + counter
  idx.PathFam.StackedHist <- c(idx.PathFam.StackedHist, list(i2))
  counter <- max(unlist(idx.PathFam.StackedHist))
}


## _ index the pathways present in each cell population----
idx.enrich.hg.StackedHist <- apply(enrich.hg, 2, function(x) which(x > 1))


## _ extract indexed frequency values from the datamatrix as a list----
enrich.hg.StackedHist <- lapply(
  idx.enrich.hg.StackedHist,
  function(x) {
    lapply(
      idx.PathFam.StackedHist,
      function(y) length(which(y %in% x)) / length(y)
    )
  }
)


## _ reformat the extracted frequency values from a list to a datamatrix----
enrich.hg.StackedHist.dm <- matrix(unlist(enrich.hg.StackedHist), ncol = ncol(enrich.hg), byrow = F)

rownames(enrich.hg.StackedHist.dm) <- names(idx.PathFam)

colnames(enrich.hg.StackedHist.dm) <- paste(colnames(enrich.hg), c("Mild", "SevCrit"), sep = "__")


## _ color palette----
col.palette <- viridis(
  nrow(enrich.hg.StackedHist.dm),
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  option = "H"
)


# _ plot Histogram to SVG file----
svg(
  "results/Stacked_Histogram.svg",
  width = 36,
  height = 14,
  pointsize = 24
)


par(fig = c(0, .8, 0, 1), mar = c(14, 6, 2, 2))
barplot(
  enrich.hg.StackedHist.dm,
  las = 2,
  col = col.palette,
  ylab = "Norm. Counts of Significantly Enriched Pathways",
  horiz = F,
  main = "Functional Enrichment Analysis",
  cex.axis = 1.25,
  cex.names = 1.5
)


par(fig = c(0.78, 1, 0, 1), mar = c(1, 1, 1, 1))
legend(
  .175,
  max(apply(enrich.hg.StackedHist.dm, 2, sum, na.rm = T)) * .9,
  legend = legend.row2,
  col = col.palette,
  pch = 15,
  pt.cex = 1.75,
  cex = .875,
  title = "Pathway Families",
  y.intersp = 1.52,
  bg = "white",
  horiz = F,
  ncol = 2,
  xpd = T
)


dev.off()


