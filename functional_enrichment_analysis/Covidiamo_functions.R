
library(heatmap3)


####_____----
## hheatmap2----
hheatmap2 <- function(pplot.hm = T, ## plot heatmap?
                      datamat0 = datamat, ## datamatrix
                      rrow.rreord = NULL, ## vector weights
                      ccol.rreord = NULL, ## vector weights
                      ddes = des, ## descriptor as data frame
                      ddes.cln = c(1, 4), ## descriptor columns to consider in the heatmap plot
                      leg.col = c(1, 4), ## column number in descriptor matrix, to be put in legend
                      xx.legend = rep(.0, 2), ## legend x position
                      yy.legend = c(.75, .928), ## legend y position
                      leg.pch = 22, ## pch legend
                      ccex.legend = .6, ## cex legend
                      ccex.pt.legend = 1, ## cex dots legend
                      llegend.y.intersp = .5, ## y interspace distance factor for legend
                      col.pal = ccol.pal, ## color palette as list
                      col.mat = colorRampPalette(c("blue", "white", "red"))(1024), ## color palette datamat
                      na20 = F, ## NA's to 0
                      mmet1 = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")[6], ## distance method
                      mmet2 = c("ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")[1], ## clustering method
                      nr.cl.row = 2, ## nr of features' clusters dimensions
                      nr.cl.col = 2, ## nr of features' clusters dimensions
                      rrow.na = .5, ## max % of NA's allowed per row in the datamatrix
                      ccol.na = .5, ## max % of NA's allowed per column in the datamatrix
                      qquan.var = .0, ## quantile threshold for selecting features
                      mmargins = c(2, 9), ## margins
                      rn.in.clmns = 1, ## rownames in columns
                      space.rep = "       ", ## space repeats to add to rownames if in columns
                      ccexRow = .5, ## cex rownames
                      ccexCol = .5, ## cex colnames
                      ddev.cur.off = F, ## dev.cur.off?
                      bbalanceColor = T, ## balanceColor?
                      sscale.hm = c("row", "column", "none")[1], ## scale
                      yes.cl.row = F, ## row color?
                      yes.cl.col = F, ## col color?
                      llog = F, ## log scale?
                      RRowv = NULL, ## vector for reordering row dendro
                      CColv = NULL, ## vector for reordering col dendro
                      sshowRowDendro = T, ## show row dendro?
                      sshowColDendro = T, ## show col dendro?
                      mmain = "", ## main title
                      blankmain = T, ## no main title
                      row.col.strip = NULL, ## row.col.strip?
                      RowLabelsColorsCln = 0, ## column of the RowColorMatrix to color row labels
                      # 	RowColorsCln=1, ## column of the RowColorMatrix to color row strip
                      nnew.plot = F, ## nnew.plot?
                      hhighlightCell = cbind(expand.grid(1:nrow(datamat.hm), 1:ncol(datamat.hm)), color = rep("transparent", nrow(datamat.hm) * ncol(datamat.hm)), lwd = rep(1, nrow(datamat.hm) * ncol(datamat.hm))), ## border color of datamat cells
                      noblankcolbars = F, ## no blank lines in colsidebars
                      noblankrowbars = F, ## no blank lines in rowsidebars
                      CColSideWidth = 1 ## height of the colsidebars
) {
  if (ddev.cur.off) dev.off(which = dev.cur())

  if (nnew.plot) x11(width = 15.2, height = 12.5, xpos = 250, ypos = 10)

  if (ncol(datamat0) > 1) datamat.hm <- data.matrix(datamat0) else datamat.hm <- data.matrix(cbind(datamat0, datamat0))

  if (na20) datamat.hm[which(is.na(datamat.hm))] <- 0

  col.na <- apply(datamat.hm, 2, function(x) length(which(is.na(x))) / length(x))
  datamat.hm <- data.matrix(datamat.hm[, which(col.na <= ccol.na)])
  row.na <- apply(datamat.hm, 1, function(x) length(which(is.na(x))) / length(x))
  datamat.hm <- data.matrix(datamat.hm[which(row.na <= rrow.na), ])

  datamat.hm.var <- apply(datamat.hm, 1, var, na.rm = T)

  wi.var <- which(datamat.hm.var >= quantile(datamat.hm.var, qquan.var, na.rm = T))
  datamat.hm <- datamat.hm[wi.var, ]

  if (!is.null(ddes)) {
    ddes.cn <- colnames(ddes)
    ddes <- as.data.frame(ddes)
    ddes <- as.data.frame(ddes[which(col.na <= ccol.na), ])
    colnames(ddes) <- ddes.cn
  }

  if (llog & min(abs(datamat.hm), na.rm = T) > 0) datamat.hm <- log(abs(datamat.hm), 10) * sign(datamat.hm)
  if (llog & min(abs(datamat.hm), na.rm = T) == 0) datamat.hm <- log(abs(datamat.hm) + min(abs(datamat.hm)[which(abs(datamat.hm) > 0)], na.rm = T) / 2, 10) * sign(datamat.hm)
  if (!llog) datamat.hm <- datamat.hm

  hcr.hm <- hclust(dist((datamat.hm), method = mmet1), method = mmet2)
  hcc.hm <- hclust(dist(t(datamat.hm), method = mmet1), method = mmet2)

  cthcr.hm <- cutree(hcr.hm, k = nr.cl.row)
  cthcc.hm <- cutree(hcc.hm, k = nr.cl.col)

  dev.set(dev.cur() - 0)

  hc.hm <- heatmap(
    datamat.hm,
    keep.dendro = T,
    col = c("blue", "white", "red"),
    scale = c("row", "column", "none")[3],
    distfun = function(x) dist(x, method = mmet1),
    hclustfun = function(x) hclust(x, method = mmet2)
  )

  ## side colors of the heatmap
  ccolnam <- ""

  if (!is.null(ddes)) {
    bblank <- as.matrix(rep(rgb(1, 1, 1), nrow(ddes)))
    ccol.hm.c0 <- as.matrix(bblank)
    ccolnam <- ""

    if (length(ddes.cln) > 0) {
      for (i1b in 1:length(ddes.cln)) {
        i1 <- ddes.cln[i1b]
        ccol.hm.c0a2 <- bblank
        ccol.hm.c0a <- ddes[, i1]

        for (i2b in 1:nlevels(as.factor(ccol.hm.c0a))) {
          i2 <- levels(as.factor(ccol.hm.c0a))[i2b]
          wi.lev <- which(ccol.hm.c0a %in% i2)
          if (length(wi.lev) > 0) ccol.hm.c0a2[wi.lev] <- col.pal[[i1b]][which(levels(as.factor(ccol.hm.c0a)) %in% i2)]
        }

        ccol.hm.c0a2 <- as.matrix(ccol.hm.c0a2)
        if (noblankcolbars) ccol.hm.c0 <- cbind(ccol.hm.c0, ccol.hm.c0a2) else ccol.hm.c0 <- cbind(ccol.hm.c0, bblank, ccol.hm.c0a2)
        ccolnam0 <- colnames(ddes)[i1]
        if (is.null(ccolnam0)) ccolnam0 <- ""
        if (noblankcolbars) ccolnam <- c(ccolnam, ccolnam0) else ccolnam <- c(ccolnam, "", ccolnam0)
      }
    }
    colnames(ccol.hm.c0) <- ccolnam
  } else {
    if (ncol(datamat0) == 1) ccol.hm.c0 <- rep(rgb(1, 1, 1), 2) else ccol.hm.c0 <- rep(rgb(1, 1, 1), ncol(datamat0))
  }

  if (!is.null(ddes)) {
    ccol.hm0b <- rep(rgb(1, 1, 1), length(cthcc.hm))
    for (i2 in unique(cthcc.hm)) ccol.hm0b[which(cthcc.hm %in% i2)] <- col.pal[[1]][i2]

    if (yes.cl.col) {
      ccol.hm0 <- cbind(ccol.hm0, bblank, ccol.hm0b)
      colnames(ccol.hm0) <- c(ccolnam, "", "clusters")
    }
  }

  rrow.hm0a <- rep(rgb(1, 1, 1), length(cthcr.hm))
  rrow.hm0 <- cbind(rrow.hm0a, rrow.hm0a)

  if (yes.cl.row) for (i2 in unique(cthcr.hm)) rrow.hm0[which(cthcr.hm %in% i2), 2] <- col.pal[[1]][i2]
  rrow.hm0 <- as.matrix(rrow.hm0)

  bblank2 <- rep(rgb(1, 1, 1), nrow(rrow.hm0))

  if (!is.null(row.col.strip)) {
    rrow.hm0 <- as.matrix(row.col.strip)
    if (!noblankrowbars & ncol(rrow.hm0) >= 2) {
      rrow.hm0b <- c()
      ccolnam.rrow.hm0 <- c()

      for (i45 in 1:ncol(rrow.hm0)) {
        rrow.hm0b <- cbind(rrow.hm0b, rrow.hm0[, i45], bblank2)
        ccolnam.rrow.hm0 <- c(ccolnam.rrow.hm0, colnames(rrow.hm0)[i45], "")
      }
    }
    colnames(rrow.hm0b) <- ccolnam.rrow.hm0

    rrow.hm0 <- rrow.hm0b
  } else {
    rrow.hm0 <- rep(rgb(1, 1, 1), length(cthcr.hm))
  }

  dev.set(dev.cur() - 0)

  if (!is.null(RRowv)) RRowv2 <- hc.hm$Rowv else RRowv2 <- NA
  if (!is.null(CColv)) CColv2 <- hc.hm$Colv else CColv2 <- NA

  datamat.hm.rn <- rownames(datamat.hm)
  if (rn.in.clmns != 0) {
    for (i in 1:min(length(datamat.hm.rn), rn.in.clmns)) {
      if (!is.null(RRowv)) rrow.ord <- hc.hm$rowInd else rrow.ord <- 1:length(datamat.hm.rn)
      sseq <- seq(i, length(datamat.hm.rn), rn.in.clmns)
      datamat.hm.rn[rrow.ord[sseq]] <- unlist(sapply(datamat.hm.rn[rrow.ord[sseq]], function(yy) paste(c(rep(space.rep, 2 * (i - 1)), yy), sep = "", collapse = "")))
    }
  } else {
    datamat.hm.rn <- rep("", nrow(datamat.hm))
  }

  rownames(datamat.hm) <- datamat.hm.rn

  if (blankmain) mmain2 <- mmain else mmain2 <- c("", "", paste(dim(datamat.hm)[], collapse = "  "), "", mmain)

  if (ncol(datamat0) == 1) colnames(datamat.hm) <- c("", "")

  ## _plot heatmap3
  hhm <- heatmap3(
    datamat.hm,
    Rowv = RRowv2,
    Colv = CColv2,
    col = col.mat,
    showRowDendro = sshowRowDendro,
    showColDendro = sshowColDendro,
    balanceColor = bbalanceColor,
    scale = sscale.hm,
    main = mmain2,
    keep.dendro = T,
    ColSideColors = ccol.hm.c0,
    ColSideWidth = CColSideWidth,
    RowSideColors = rrow.hm0,
    RowAxisColors = RowLabelsColorsCln,
    margins = mmargins,
    cexRow = ccexRow,
    cexCol = ccexCol,
    na.rm = T,
    highlightCell = hhighlightCell
  )

  ## _legend
  if (!is.null(leg.col) & !is.null(ddes)) {
    for (i7 in 1:length(leg.col)) {
      legend(
        x = xx.legend[i7],
        y = yy.legend[i7],
        legend = levels(as.factor(ddes[, leg.col[i7]])),
        pch = leg.pch,
        pt.cex = ccex.pt.legend,
        col = col.pal[[i7]][1:nlevels(as.factor(ddes[, leg.col[i7]]))],
        cex = ccex.legend,
        title = colnames(ddes)[leg.col[i7]],
        y.intersp = llegend.y.intersp,
        title.adj = 0
      )
    }
  }

  return(hhm)
}

	

	
	


