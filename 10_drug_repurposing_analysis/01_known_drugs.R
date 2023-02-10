library(pheatmap)
setwd("working directory")
source("covidiamo_functions_RAGE.R")



data_trt<-get_data(trt="trt_cp")
gene_ind<-get_gene_ind(all_data=data_trt)
enrichment<-enrich_data(all_data=data_trt, ind=gene_ind)
ids<-enrich_trt(all_enrich=enrichment, trt="trt_cp", to_plot=T)

#Bcells
enrichment[[1]][[1]][c("dexamethasone", "baricitinib", "ritonavir")]
#T cells
enrichment[[2]][[1]][c("dexamethasone", "baricitinib", "ritonavir")]
#NK cells
enrichment[[3]][[1]][c("dexamethasone", "baricitinib", "ritonavir")]
#Myeloid cells
enrichment[[4]][[1]][c("dexamethasone", "baricitinib", "ritonavir", "azithromycin", "hydroxychloroquine")]
#CP
enrichment[[5]][[1]][c("dexamethasone", "baricitinib", "ritonavir", "azithromycin", "hydroxychloroquine")]


#########
##Randomization baricitinib
#########

set.seed(194876)
ind<-gene_ind[[1]]
cMAP<-data_trt[[1]]
Bcell_rand<-list()
for(iter in 1:10000){
ind_rand<-sample(1:nrow(cMAP), length(ind), replace = F)
if(sum(cMAP$pert_iname[ind_rand]=="baricitinib")>1){
  a<-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
  c<-sum(cMAP$pert_iname=="baricitinib")-a
  b<-length(ind_rand)-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
  d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname=="baricitinib"),ind_rand )))
  Bcell_rand[[iter]]<-fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]]
}
}

sum(unlist(Bcell_rand)<enrichment[[1]][[1]]["baricitinib"])

set.seed(194876)
ind<-gene_ind[[3]]
cMAP<-data_trt[[2]]
Tcell_rand<-list()
for(iter in 1:10000){
  ind_rand<-sample(1:nrow(cMAP), length(ind), replace = F)
  if(sum(cMAP$pert_iname[ind_rand]=="baricitinib")>1){
    a<-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    c<-sum(cMAP$pert_iname=="baricitinib")-a
    b<-length(ind_rand)-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname=="baricitinib"),ind_rand )))
    Tcell_rand[[iter]]<-fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]]
  }
}

sum(unlist(Tcell_rand)<enrichment[[2]][[1]]["baricitinib"])


set.seed(194876)
ind<-gene_ind[[5]]
cMAP<-data_trt[[3]]
NKell_rand<-list()
for(iter in 1:10000){
  ind_rand<-sample(1:nrow(cMAP), length(ind), replace = F)
  if(sum(cMAP$pert_iname[ind_rand]=="baricitinib")>1){
    a<-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    c<-sum(cMAP$pert_iname=="baricitinib")-a
    b<-length(ind_rand)-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname=="baricitinib"),ind_rand )))
    NKell_rand[[iter]]<-fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]]
  }
}

sum(unlist(NKell_rand)<enrichment[[3]][[1]]["baricitinib"])


set.seed(194876)
ind<-gene_ind[[7]]
cMAP<-data_trt[[4]]
Myeloid_rand<-list()
for(iter in 1:10000){
  ind_rand<-sample(1:nrow(cMAP), length(ind), replace = F)
  if(sum(cMAP$pert_iname[ind_rand]=="baricitinib")>1){
    a<-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    c<-sum(cMAP$pert_iname=="baricitinib")-a
    b<-length(ind_rand)-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname=="baricitinib"),ind_rand )))
    Myeloid_rand[[iter]]<-fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]]
  }
}

sum(unlist(Myeloid_rand)<enrichment[[4]][[1]]["baricitinib"])


set.seed(194876)
ind<-gene_ind[[9]]
cMAP<-data_trt[[5]]
CP_rand<-list()
for(iter in 1:10000){
  ind_rand<-sample(1:nrow(cMAP), length(ind), replace = F)
  if(sum(cMAP$pert_iname[ind_rand]=="baricitinib")>1){
    a<-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    c<-sum(cMAP$pert_iname=="baricitinib")-a
    b<-length(ind_rand)-sum(cMAP$pert_iname[ind_rand]=="baricitinib")
    d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname=="baricitinib"),ind_rand )))
    CP_rand[[iter]]<-fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]]
  }
}

sum(unlist(CP_rand)<enrichment[[5]][[1]]["baricitinib"])