library(slinky)

# update following lines with your details:
user_key <- "personal user key"
gctx <- "path to level4 data"
info <- "path to metadata"
sl <- Slinky(user_key, gctx, info)


##load info about perturbagens
sig<-read.csv(info, sep="\t")
col.ix <- which(sig$pert_iname == "baricitinib")


data <- readGCTX(sl[, col.ix-129543])#match col indexes 
gene_info<-read.csv(info, sep="\t")
data_mean<-rowMeans(data)

library(openxlsx)
IPA<-read.xlsx("data/PredictionIpa.xlsx",1)
RAGE_path<-IPA[,1]

names(data_mean)<-gene_info$pr_gene_symbol[match(names(data_mean), gene_info$pr_gene_id)]

write.csv(cbind(names(data_mean[RAGE_path]), data_mean[RAGE_path]), file="baricitinib_zscore.csv", row.names = T)
deact<-IPA[IPA$Expr.Log.Ratio<0,1]
act<-IPA[IPA$Expr.Log.Ratio>0,1]

wilcox.test(data_mean[deact], data_mean, alternative="less")

