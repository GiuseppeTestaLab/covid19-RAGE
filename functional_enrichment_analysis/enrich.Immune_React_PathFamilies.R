
setwd("functional_enrichment_analysis")
source("Covidiamo_functions.R")

####_____----
## upload and filter functional enrichment results----

# _ upload functional enrichment results----
funct.enrich <- read.table("data/funct.enrich.tsv", sep = "\t", header = T, quote='')


# _ select entries in funct.enrich by FDR, max number of background genes and min number of genes per enriched functional term----
idx.enrich.1a <- which((
  as.numeric(as.vector(funct.enrich$fdr)) <= 0.001 &
    as.numeric(as.vector(funct.enrich$number_of_genes_in_background)) <= 500 &
    as.numeric(as.vector(funct.enrich$number_of_genes)) >= 3
))


# _ select entries in funct.enrich by dual trend patterns----
idx.enrich.dual_patterns <- unique(unlist(sapply(
  c("1\\.1", "1\\.2", "1\\.3", "2\\.1", "2\\.2", "2\\.3", "3\\.1", "3\\.2", "3\\.3"),
  function(x) {
    grep(x, funct.enrich$ID)
  }
)))


# _ select entries in funct.enrich by single trend patterns----
idx.enrich.single_patterns <- setdiff(1:nrow(funct.enrich), idx.enrich.dual_patterns)


# _ select entries in funct.enrich by "Mild" or "SevCrit" severity----
idx.enrich.2a <- unique(unlist(
  sapply(
    c("Mild", "SevCrit"),
    function(x) {
      grep(x, funct.enrich[, 1])
    }
  )
))


# _ select entries in funct.enrich by increasing or decreasing trend of expression----
idx.enrich.2b <- unique(unlist(
  sapply(
    c("123", "321"),
    function(x) {
      grep(x, funct.enrich[, 1])
    }
  )
))


# _ select entries in funct.enrich by database----
idx.enrich.2c <- grep("hsa", funct.enrich[, 3], invert = T)


# _ select entries in funct.enrich intersecting all the above criteria----
idx.enrich <- Reduce(
  "intersect",
  lapply(
    list(
      idx.enrich.1a,
      idx.enrich.single_patterns,
      idx.enrich.2a,
      idx.enrich.2b,
      idx.enrich.2c
    ),
    function(x) x
  )
)


# _ merge GO IDs and GO description----
termID <- apply(funct.enrich[, c(3, 11)], 1, paste, sep = ": ", collapse = ": ")


# _ filtered funct.enrich----
funct.enrich.02 <- funct.enrich[idx.enrich, ]


# _ re-format entries (for better readability) in some fields of filtered funct.enrich.02----
funct.enrich.02$ID <- gsub("\\.", "_", funct.enrich.02$ID)

funct.enrich.02$number_of_genes <- gsub(" ", "", funct.enrich.02$number_of_genes)

funct.enrich.02$number_of_genes_in_background <- gsub(" ", "", funct.enrich.02$number_of_genes_in_background)



####_____----
## prepare cell families, patterns and metadata info----

# _ cell families and patterns----
funct.enrich.02.famnam <- as.data.frame(t(sapply(
  funct.enrich.02$ID,
  function(x) strsplit(x, "__")[[1]]
)))

colnames(funct.enrich.02.famnam) <- c("CellFam", "Pattern")


# _ merge GO IDs and GO description in the filtered funct.enrich.02 ----
termID.02 <- apply(funct.enrich.02[, c(3, 11)], 1, paste, sep = ": ", collapse = ": ")


# _ frequencies table of GO terms and cells by patterns----
funct.enrich.02.ta <- table(
  funct.enrich.02$ID,
  termID.02
)


# _ set up a data matrix w/ logged FDR values to be visualized as HeatMap----
funct.enrich.02.ta.pv <- funct.enrich.02.ta

funct.enrich.02.ta.pv[which(funct.enrich.02.ta.pv >= 0)] <- 1


# _ populating the data matrix w/ logged FDR values----
for (i in rownames(funct.enrich.02.ta)) {
  wix1 <- which(funct.enrich.02$ID %in% i)
  wix2r <- which(rownames(funct.enrich.02.ta) %in% i)
  x1 <- names(which(funct.enrich.02.ta[wix2r, ] > 0))
  wix2c <- which(colnames(funct.enrich.02.ta) %in% x1)
  wix3 <- which(termID.02[wix1] %in% x1)
  for (i2 in wix2c) {
    wix4 <- which(termID.02[wix1[wix3]] %in% colnames(funct.enrich.02.ta)[i2])
    funct.enrich.02.ta.pv[wix2r, i2] <- log(-log(min(as.numeric(funct.enrich.02$fdr[wix1[wix3][wix4]]), na.rm = T), 10), 2)
  }
}


# _ order rows and columns of frequencies table of GO terms and cells by patterns----
funct.enrich.02.ta.rn2 <- sapply(
  rownames(funct.enrich.02.ta),
  function(x) strsplit(x, "__")[[1]][1]
)


idx.famnam <- unlist(sapply(
  unique(funct.enrich.02.famnam$CellFam),
  function(x) which(funct.enrich.02.ta.rn2 %in% x)
))


funct.enrich.02.ta2 <- funct.enrich.02.ta.pv[idx.famnam, ]

wi.app.1 <- which(apply(funct.enrich.02.ta2, 1, sum, na.rm = T) > 0)

wi.app.2 <- which(apply(funct.enrich.02.ta2, 2, sum, na.rm = T) > 0)


funct.enrich.02.ta3 <- funct.enrich.02.ta2[wi.app.1, wi.app.2]


famnam2 <- sapply(
  rownames(funct.enrich.02.ta3),
  function(x) strsplit(x, "__")[[1]][2]
)

rownames(funct.enrich.02.ta3) <- names(famnam2)



####_____----
## setup data and metadata matrices----

# _ setup and format the data matrix w/ logged FDR values for HeatMap----
datamat.hm <- as.matrix(t(funct.enrich.02.ta3))


# _ prepare metadata matrix for HeatMap----
metadatamat.hm <- DataFrame(t(as.matrix(
  sapply(
    colnames(datamat.hm),
    function(x) strsplit(x, "__")[[1]]
  )
)))

colnames(metadatamat.hm) <- c("CellFam.0", "CellFamSev")


# _ split first field of metadatamat.hm into severity and cells info----
ds.tmp1 <- sapply(metadatamat.hm$CellFam.0, function(x) strsplit(x, "_")[[1]][1])

ds.tmp2 <- sapply(metadatamat.hm$CellFam.0, function(x) strsplit(x, "_")[[1]][2])


# _ add severity and cells population info into metadatamat.hm----
metadatamat.hm <- cbind(
  metadatamat.hm,
  Severity = ds.tmp1,
  CellFam. = ds.tmp2
)

## _ reorder datamat.hm0 and metadatamat.hm entries by severity----
oord.hm <- order(metadatamat.hm$Severity)

metadatamat.hm <- metadatamat.hm[oord.hm, ]

datamat.hm <- datamat.hm[, oord.hm]


## _ group GO terms in bigger Immune-Response GO terms families----
RegImmResp <- c("GO.0006954", "GO.0045088", "GO.0035456", "GO.0002697", "GO.0034340", "GO.0046598", "GO.0050688", "GO.0032479", "GO.0002825", "GO.0001818", "GO.0032480", "GO.0046596", "GO.0050727", "GO.0002698", "GO.0042742", "GO.0050777", "GO.0006959", "GO.0045089", "GO.0045824", "GO.0002683", "GO.0050786")
AntigProc <- c("GO.0042605", "GO.0019886", "GO.0019884", "GO.0032395", "GO.0042613")
IFNPathAndViralResp <- c("GO.0060333", "GO.0050792", "GO.0048525", "GO.1903900", "GO.0043903", "GO.0045069", "GO.1903901", "GO.0035455", "GO.0034341", "GO.0071346", "GO.0051607", "GO.0009615", "GO.0045071", "GO.0060337")


BigPathFam <- list(
  RegImmResp,
  AntigProc,
  IFNPathAndViralResp
)

names(BigPathFam) <- c("Regulation of Immune Response", "Antigen Processing", "IFN Pathway and Viral Response")


# _ index the GO terms associated to the bigger GO terms families within datamat ----
idx.BigPathFam <- lapply(
  BigPathFam,
  function(x) {
    unlist(sapply(
      x,
      function(y) grep(y, rownames(datamat.hm))
    ))
  }
)



####_____----
## prepare data and metadata matrices for plotting----

# _ define colors for heatmap metadata and legends----
hm.row.col <- rep(rgb(1, 1, 1), nrow(datamat.hm))

hm.row.col <- cbind(hm.row.col, hm.row.col, hm.row.col)

legend.row <- names(idx.BigPathFam)

legend.row.col <- rep(rgb(1, 1, 1), length(idx.BigPathFam))


for (i in 1:length(idx.BigPathFam)) {
  legend.row.col[i] <-
    hm.row.col[idx.BigPathFam[[i]], 1] <-
    rgb(i / (length(idx.BigPathFam) + 1), i / (length(idx.BigPathFam) + 1), 1)
}


# _ subset data and metadata colors matrices by bigger Immune-Response GO terms families----
datamat.hm.02 <- datamat.hm[unlist(idx.BigPathFam), ]

hm.row.col.02 <- hm.row.col[unlist(idx.BigPathFam), ]


# _ reformat colnames of datamatrix for better readability----
datamat.hm.02.cn <- colnames(datamat.hm.02)
datamat.hm.02.cn <- gsub("__", "@", datamat.hm.02.cn)
datamat.hm.02.cn <- gsub("Myeloid", "Myel.", datamat.hm.02.cn)
datamat.hm.02.cn <- gsub("Mild_", "", datamat.hm.02.cn)
datamat.hm.02.cn <- gsub("SevCrit_", "", datamat.hm.02.cn)
datamat.hm.02.cn <- gsub("@", "_", datamat.hm.02.cn)

gsubin <- c("1_1", "1_2", "1_3", "2_2", "2_3", "3_3", "2_1", "3_1", "3_2")
gsubout <- c("1.1", "1.2", "1.3", "2.2", "2.3", "3.3", "2.1", "3.1", "3.2")

for (i in 1:length(gsubin)) datamat.hm.02.cn <- gsub(gsubin[i], gsubout[i], datamat.hm.02.cn)

colnames(datamat.hm.02) <- datamat.hm.02.cn


## _ reorder cells families in datamat and metadatamat----
ord.hm.cn <- order(
  metadatamat.hm$CellFam.,
  colnames(datamat.hm.02)
)

datamat.hm.02 <- datamat.hm.02[, ord.hm.cn]

metadatamat.hm.02 <- metadatamat.hm[ord.hm.cn, ]


# _ reformat data and metadata matrices rows and columns for plotting----
datamat.hm.02.nr <- nrow(datamat.hm.02)

datamat.hm.02.cn <- colnames(datamat.hm.02)


datamat.hm.02 <- cbind(
  datamat.hm.02[, 1:3],
  rep(1, datamat.hm.02.nr),
  datamat.hm.02[, 4],
  rep(1, datamat.hm.02.nr),
  datamat.hm.02[, 5:8],
  datamat.hm.02[, 9:14]
)


datamat.hm.02.cn <- c(
  datamat.hm.02.cn[1:3],
  datamat.hm.02.cn[3],
  datamat.hm.02.cn[4],
  datamat.hm.02.cn[4],
  datamat.hm.02.cn[5:8],
  datamat.hm.02.cn[9:14]
)

colnames(datamat.hm.02) <- datamat.hm.02.cn


metadatamat.hm.02 <- rbind(
  metadatamat.hm.02[1:3, ],
  metadatamat.hm.02[3, ],
  metadatamat.hm.02[4, ],
  metadatamat.hm.02[4, ],
  metadatamat.hm.02[5:8, ],
  metadatamat.hm.02[9:14, ]
)

metadatamat.hm.02$Severity[c(4, 6)] <- c("SevCrit", "SevCrit")


idx.NonEmpty <- Reduce("intersect", lapply(
  sapply(
    c("CP_123", "NK_123"),
    function(y) list(grep(y, colnames(datamat.hm.02), invert = T))
  ),
  function(x) x
))


datamat.hm.02 <- datamat.hm.02[, idx.NonEmpty]

metadatamat.hm.02 <- metadatamat.hm.02[idx.NonEmpty, ]

colnames(metadatamat.hm.02)[3] <- ""


# _ scale datamat----
datamat.hm.03 <- datamat.hm.02 / max(datamat.hm.02)


# _ count nr. of genes in each functional enrichment term associated w/ the RAGE pathway----
RagePathGenes <- unlist(lapply(
  strsplit(funct.enrich$preferredNames, ","),
  function(x) length(which(x %in% c("RAGE", "AGER", "FPR1", "HMGB1", "HMGB2", "S100A12", "S100A13", "S100A4", "S100A7", "S100A8", "S100A9", "S100B")))
))


# _ select functional enrichment terms w/ at least 1 gene associated w/ the RAGE pathway----
idx.RagePathGenes <- which(RagePathGenes >= 1)

idx.RagePathGenes <- intersect(idx.enrich, idx.RagePathGenes)


# _ index functional enrichment terms w/ associated to RAGE pathway within the datamat----
idx.rage.hm <- which(rownames(datamat.hm.03) %in% termID[idx.RagePathGenes])


# _ define RAGE pathways colors----
hm.row.col.02[idx.rage.hm, 3] <- rgb(1, .65, 0)

hm.row.col.02[grep("RAGE", rownames(datamat.hm.03)), 3] <- rgb(1, .25, 0)

colnames(hm.row.col.02) <- rep("", ncol(hm.row.col.02))

RRowv0 <- x.ord2 <- sapply(hm.row.col.02[, 1], function(x) which(unique(hm.row.col.02[, 1]) %in% x))


# _ plot Heatmap to SVG file----
svg(
  "results/Heatmap.svg",
  width = 28,
  height = 16,
  pointsize = 24
)

par(fig = c(0, 1, 0, 1), new = F, mar = c(0, 0, 0, 0), bg = "white")

hheatmap.out2 <- hheatmap2(
  pplot.hm = T,
  datamat0 = datamat.hm.03,
  rrow.rreord = x.ord2,
  ccol.rreord = NULL,
  ddes = metadatamat.hm.02,
  ddes.cln = 3,
  leg.col = 3,
  xx.legend = 1.815,
  yy.legend = 1.825,
  ccex.legend = .75,
  ccex.pt.legend = 1.75,
  llegend.y.intersp = 1.2,
  col.pal = list(c(rgb(0.7, 0.7, 0.7), rgb(0, 0, 0))),
  col.mat = colorRampPalette(
    c(
      rep(rgb(1, 1, 1), 2),
      rep(rgb(1, 0, 0), 1)
    ),
    bias = 1.12,
    interpolate = "linear"
  )(1024),
  na20 = F,
  mmet1 = "manhattan",
  mmet2 = "ward.D2",
  nr.cl.row = 2,
  nr.cl.col = 2,
  rrow.na = .5,
  ccol.na = .5,
  qquan.var = .0,
  mmargins = c(8, 16),
  rn.in.clmns = 1,
  ccexRow = 0.8,
  ccexCol = 1.6,
  ddev.cur.off = F,
  bbalanceColor = T,
  sscale.hm = "none",
  yes.cl.row = F,
  yes.cl.col = F,
  llog = F,
  RRowv = RRowv0,
  CColv = NULL,
  sshowRowDendro = F,
  sshowColDendro = F,
  mmain = "",
  blankmain = T,
  row.col.strip = hm.row.col.02,
  RowLabelsColorsCln = 0,
  nnew.plot = F,
  hhighlightCell = cbind(expand.grid(1:nrow(datamat.hm.03), 1:ncol(datamat.hm.03)), color = rep(rgb(0, 0, 0, .085), nrow(datamat.hm.03) * ncol(datamat.hm.03)), lwd = rep(1.0, nrow(datamat.hm.03) * ncol(datamat.hm.03))),
  noblankcolbars = F,
  noblankrowbars = F,
  CColSideWidth = 1
)


legend(
  .05, .45,
  legend = (names(BigPathFam)),
  col = (unique(hm.row.col.02[, 1])),
  pch = 15,
  pt.cex = 1.75,
  cex = .875,
  title = NULL,
  y.intersp = 1.2,
  bg = "white"
)


legend(
  .05, .25,
  legend = c(paste("At least 1 RAGE Binding Gene", sep = "", collapse = ""), "RAGE Binding Genes (GO:0050786)"),
  col = c(rgb(1, .65, 0), rgb(1, .25, 0)),
  pch = 15,
  pt.cex = 1.75,
  cex = .875,
  title = NULL,
  y.intersp = 1.2,
  bg = "white",
  trace = T
)


legend(
  .05, .65,
  legend = c("Mild", "SevCrit"),
  col = c(rgb(0.7, 0.7, 0.7), rgb(0, 0, 0)),
  pch = 15,
  pt.cex = 1.75,
  cex = .875,
  title = "Severity",
  y.intersp = 1.2,
  bg = "white"
)


text(.27, .817, sscale.legend, font = 1, cex = .95)


dev.off()




