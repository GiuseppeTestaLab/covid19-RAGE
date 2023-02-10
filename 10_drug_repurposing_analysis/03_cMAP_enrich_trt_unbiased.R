#################
### Enrichment
################
library(pheatmap)

enrich<-function(cMAP=Bcell, ind=ind_up){
drug_Bcell_up<-names(which(table(cMAP$pert_iname[ind])>1))
drug_Bcell_up_p<-c()
for(drug in drug_Bcell_up){
  a<-sum(cMAP$pert_iname[ind]==drug)
  c<-sum(cMAP$pert_iname==drug)-a
  b<-length(ind)-sum(cMAP$pert_iname[ind]==drug)
  d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname==drug),ind )))
  drug_Bcell_up_p<-c(drug_Bcell_up_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
}
names(drug_Bcell_up_p)<-drug_Bcell_up

moa_ind<-unlist(strsplit(cMAP$moa[ind], "\\|"))
moa_all<-unlist(strsplit(cMAP$moa, "\\|"))
moa_Bcell_up<-names(which(table(moa_ind)>1))
moa_Bcell_up_p<-c()
for(moa in moa_Bcell_up){
  a<-sum(moa_ind==moa)
  c<-sum(moa_all==moa)-a
  b<-length(moa_ind)-sum(moa_ind==moa)
  moa_back<-unlist(strsplit(cMAP$moa[-ind], "\\|"))
  d<-length(moa_back[moa_back!=moa])
  moa_Bcell_up_p<-c(moa_Bcell_up_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
}
names(moa_Bcell_up_p)<-moa_Bcell_up


target_ind<-unlist(strsplit(cMAP$target_name[ind], "\\|"))
target_all<-unlist(strsplit(cMAP$target_name, "\\|"))
target_Bcell_up<-names(which(table(target_ind)>1))
target_Bcell_up_p<-c()
for(target in target_Bcell_up){
  a<-sum(target_ind==target)
  c<-sum(target_all==target)-a
  b<-length(target_ind)-sum(target_ind==target)
  target_back<-unlist(strsplit(cMAP$target_name[-ind], "\\|"))
  d<-length(moa_back[moa_back!=target])
  target_Bcell_up_p<-c(target_Bcell_up_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
}
names(target_Bcell_up_p)<-target_Bcell_up

pvals<-list(drug=drug_Bcell_up_p, moa=moa_Bcell_up_p, target=target_Bcell_up_p)
return(pvals)
}



get_data<-function(trt){
Bcell<-read.csv("data/Bcell.gct", sep="\t", header=T)
Bcell<-Bcell[-1,]
Bcell<-Bcell[Bcell$pert_type==trt,]


Tcell<-read.csv("data/Tcell.gct",  sep="\t", header=T)
Tcell<-Tcell[-1,]
Tcell<-Tcell[Tcell$pert_type==trt,]


NKcell<-read.csv("data/NKcell.gct", sep="\t", header=T)
NKcell<-NKcell[-1,]
NKcell<-NKcell[NKcell$pert_type==trt,]


Myeloid<-read.csv("data/Myeloid.gct", sep="\t", header=T)
Myeloid<-Myeloid[-1,]
Myeloid<-Myeloid[Myeloid$pert_type==trt,]

all_data<-list(Bcell, Tcell, NKcell, Myeloid)
return(all_data)
}


enrich_data<-function(all_data){
  Bcell<-all_data[[1]]
  ind_up<-which(Bcell$fdr_q_nlog10>2 & Bcell$raw_cs>0)
  ind_dn<-which(Bcell$fdr_q_nlog10>2 & Bcell$raw_cs<0)
  Bcell_p_up<-enrich(cMAP=Bcell, ind=ind_up)
  Bcell_p_dn<-enrich(cMAP=Bcell, ind=ind_dn)
  
  
  Tcell<-all_data[[2]]
  ind_up<-which(Tcell$fdr_q_nlog10>2 & Tcell$raw_cs>0)
  ind_dn<-which(Tcell$fdr_q_nlog10>2 & Tcell$raw_cs<0)
  Tcell_p_up<-enrich(cMAP=Tcell, ind=ind_up)
  Tcell_p_dn<-enrich(cMAP=Tcell, ind=ind_dn)
  
  NKcell<-all_data[[3]]
  ind_up<-which(NKcell$fdr_q_nlog10>2 & NKcell$raw_cs>0)
  ind_dn<-which(NKcell$fdr_q_nlog10>2 & NKcell$raw_cs<0)
  NKcell_p_up<-enrich(cMAP=NKcell, ind=ind_up)
  NKcell_p_dn<-enrich(cMAP=NKcell, ind=ind_dn)
  
  
  Myeloid<-all_data[[4]]
  ind_up<-which(Myeloid$fdr_q_nlog10>2 & Myeloid$raw_cs>0)
  ind_dn<-which(Myeloid$fdr_q_nlog10>2 & Myeloid$raw_cs<0)
  Myeloid_p_up<-enrich(cMAP=Myeloid, ind=ind_up)
  Myeloid_p_dn<-enrich(cMAP=Myeloid, ind=ind_dn)
  all_enrich<-list(Bcell_p_up, Tcell_p_up, NKcell_p_up, Myeloid_p_up,
                 Bcell_p_dn, Tcell_p_dn, NKcell_p_dn, Myeloid_p_dn)
  return(all_enrich)
}


enrich_trt<-function(all_enrich, trt, thr_up=2,thr_up_target=1, thr_dn=1){
  Bcell_p_up<-all_enrich[[1]]
  Tcell_p_up<-all_enrich[[2]]
  NKcell_p_up<-all_enrich[[3]]
  Myeloid_p_up<-all_enrich[[4]]
  Bcell_p_dn<-all_enrich[[5]]
  Tcell_p_dn<-all_enrich[[6]]
  NKcell_p_dn<-all_enrich[[7]]
  Myeloid_p_dn<-all_enrich[[8]]
  
  drugs_up<-c()
  targets_up<-c()
  moas_up<-c()
  drugs_dn<-c()
  targets_dn<-c()
  moas_dn<-c()
  
drugs_Bcell_sel<-names(Bcell_p_up[[1]])[p.adjust(Bcell_p_up[[1]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_up[[1]])[p.adjust(Tcell_p_up[[1]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_up[[1]])[p.adjust(NKcell_p_up[[1]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_up[[1]])[p.adjust(Myeloid_p_up[[1]], method="fdr")<=0.05]
drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=4,
                  dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_up[[1]], method="fdr")[drugs_allcells_up]
drug_heat[,2]<-p.adjust(Tcell_p_up[[1]], method="fdr")[drugs_allcells_up]
drug_heat[,3]<-p.adjust(NKcell_p_up[[1]], method="fdr")[drugs_allcells_up]
drug_heat[,4]<-p.adjust(Myeloid_p_up[[1]], method="fdr")[drugs_allcells_up]

drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>thr_up,]

drug_heat<-(-log10(drug_heat))


if(!is.null(nrow(drug_heat))){
  drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
  drugs_up<-rownames(drug_heat)
  
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(drug_heat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_drugs_up",trt,"enrich.png",sep="_"), res=600, 2000, 4000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}

if(trt %in% c("trt_cp", "ATC")){

#################
##moa
###################

drugs_Bcell_sel<-names(Bcell_p_up[[2]])[p.adjust(Bcell_p_up[[2]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_up[[2]])[p.adjust(Tcell_p_up[[2]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_up[[2]])[p.adjust(NKcell_p_up[[2]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_up[[2]])[p.adjust(Myeloid_p_up[[2]], method="fdr")<=0.05]
drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=4,
                  dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_up[[2]], method="fdr")[drugs_allcells_up]
drug_heat[,2]<-p.adjust(Tcell_p_up[[2]], method="fdr")[drugs_allcells_up]
drug_heat[,3]<-p.adjust(NKcell_p_up[[2]], method="fdr")[drugs_allcells_up]
drug_heat[,4]<-p.adjust(Myeloid_p_up[[2]], method="fdr")[drugs_allcells_up]

drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>1,]
drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
drug_heat<-(-log10(drug_heat))

if(!is.null(nrow(drug_heat))){
drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]

moas_up<-rownames(drug_heat)


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_moa_up",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}

#################
##targets
###################

drugs_Bcell_sel<-names(Bcell_p_up[[3]])[p.adjust(Bcell_p_up[[3]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_up[[3]])[p.adjust(Tcell_p_up[[3]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_up[[3]])[p.adjust(NKcell_p_up[[3]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_up[[3]])[p.adjust(Myeloid_p_up[[3]], method="fdr")<=0.05]
drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=4,
                  dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_up[[3]], method="fdr")[drugs_allcells_up]
drug_heat[,2]<-p.adjust(Tcell_p_up[[3]], method="fdr")[drugs_allcells_up]
drug_heat[,3]<-p.adjust(NKcell_p_up[[3]], method="fdr")[drugs_allcells_up]
drug_heat[,4]<-p.adjust(Myeloid_p_up[[3]], method="fdr")[drugs_allcells_up]

drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>thr_up_target,]
drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
drug_heat<-(-log10(drug_heat))
if(!is.null(nrow(drug_heat))){
  
drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]

targets_up<-rownames(drug_heat)


paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_target_up",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}

if(trt=="ATC"){
  #################
  ##targets
  ###################
  
  drugs_Bcell_sel<-names(Bcell_p_up[[4]])[p.adjust(Bcell_p_up[[4]], method="fdr")<=0.05]
  drugs_Tcell_sel<-names(Tcell_p_up[[4]])[p.adjust(Tcell_p_up[[4]], method="fdr")<=0.05]
  drugs_NKcell_sel<-names(NKcell_p_up[[4]])[p.adjust(NKcell_p_up[[4]], method="fdr")<=0.05]
  drugs_Myeloid_sel<-names(Myeloid_p_up[[4]])[p.adjust(Myeloid_p_up[[4]], method="fdr")<=0.05]
  drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))
  
  drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=4,
                    dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid")))
  
  drug_heat[,1]<-p.adjust(Bcell_p_up[[4]], method="fdr")[drugs_allcells_up]
  drug_heat[,2]<-p.adjust(Tcell_p_up[[4]], method="fdr")[drugs_allcells_up]
  drug_heat[,3]<-p.adjust(NKcell_p_up[[4]], method="fdr")[drugs_allcells_up]
  drug_heat[,4]<-p.adjust(Myeloid_p_up[[4]], method="fdr")[drugs_allcells_up]
  
  drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>thr_up_target,]
  drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
  drug_heat<-(-log10(drug_heat))
  if(!is.null(nrow(drug_heat))){
    
    drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
    
    targets_up<-rownames(drug_heat)
    
    
    paletteLength <- 50
    # use floor and ceiling to deal with even/odd length pallettelengths
    myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
    length(myBreaks) == length(paletteLength) + 1
    
    png(paste("results/heat_fourth_up",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
    print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
    dev.off()
  }
  
}
}
############################################DN###############################
#################
##drugs
###################

drugs_Bcell_sel<-names(Bcell_p_dn[[1]])[p.adjust(Bcell_p_dn[[1]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_dn[[1]])[p.adjust(Tcell_p_dn[[1]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_dn[[1]])[p.adjust(NKcell_p_dn[[1]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_dn[[1]])[p.adjust(Myeloid_p_dn[[1]], method="fdr")<=0.05]
drugs_allcells_dn<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_dn), ncol=4,
                  dimnames=list(c(drugs_allcells_dn), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_dn[[1]], method="fdr")[drugs_allcells_dn]
drug_heat[,2]<-p.adjust(Tcell_p_dn[[1]], method="fdr")[drugs_allcells_dn]
drug_heat[,3]<-p.adjust(NKcell_p_dn[[1]], method="fdr")[drugs_allcells_dn]
drug_heat[,4]<-p.adjust(Myeloid_p_dn[[1]], method="fdr")[drugs_allcells_dn]

drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>thr_dn,]

drug_heat<-(-log10(drug_heat))

if(!is.null(nrow(drug_heat)) & nrow(drug_heat)>0){
  drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
  drugs_dn<-rownames(drug_heat)
  
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(drug_heat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_drugs_dn",trt,"enrich.png",sep="_"), res=600, 2000, 4000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}

if(trt %in% c("trt_cp", "ATC")){
#################
##moa
###################

drugs_Bcell_sel<-names(Bcell_p_dn[[2]])[p.adjust(Bcell_p_dn[[2]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_dn[[2]])[p.adjust(Tcell_p_dn[[2]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_dn[[2]])[p.adjust(NKcell_p_dn[[2]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_dn[[2]])[p.adjust(Myeloid_p_dn[[2]], method="fdr")<=0.05]
drugs_allcells_dn<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_dn), ncol=4,
                  dimnames=list(c(drugs_allcells_dn), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_dn[[2]], method="fdr")[drugs_allcells_dn]
drug_heat[,2]<-p.adjust(Tcell_p_dn[[2]], method="fdr")[drugs_allcells_dn]
drug_heat[,3]<-p.adjust(NKcell_p_dn[[2]], method="fdr")[drugs_allcells_dn]
drug_heat[,4]<-p.adjust(Myeloid_p_dn[[2]], method="fdr")[drugs_allcells_dn]

#drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>2,]
drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
drug_heat<-(-log10(drug_heat))

if(!is.null(nrow(drug_heat))){
  drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
  
  moas_dn<-rownames(drug_heat)
  
  
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_moa_dn",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}


#################
##targets
###################

drugs_Bcell_sel<-names(Bcell_p_dn[[3]])[p.adjust(Bcell_p_dn[[3]], method="fdr")<=0.05]
drugs_Tcell_sel<-names(Tcell_p_dn[[3]])[p.adjust(Tcell_p_dn[[3]], method="fdr")<=0.05]
drugs_NKcell_sel<-names(NKcell_p_dn[[3]])[p.adjust(NKcell_p_dn[[3]], method="fdr")<=0.05]
drugs_Myeloid_sel<-names(Myeloid_p_dn[[3]])[p.adjust(Myeloid_p_dn[[3]], method="fdr")<=0.05]
drugs_allcells_dn<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))

drug_heat<-matrix(nrow=length(drugs_allcells_dn), ncol=4,
                  dimnames=list(c(drugs_allcells_dn), c("Bcell", "Tcell", "NKcell", "Myeloid")))

drug_heat[,1]<-p.adjust(Bcell_p_dn[[3]], method="fdr")[drugs_allcells_dn]
drug_heat[,2]<-p.adjust(Tcell_p_dn[[3]], method="fdr")[drugs_allcells_dn]
drug_heat[,3]<-p.adjust(NKcell_p_dn[[3]], method="fdr")[drugs_allcells_dn]
drug_heat[,4]<-p.adjust(Myeloid_p_dn[[3]], method="fdr")[drugs_allcells_dn]

#drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>2,]
drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
drug_heat<-(-log10(drug_heat))

if(nrow(drug_heat)>0){
  drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
  
  targets_dn<-rownames(drug_heat)
  
  

paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
length(myBreaks) == length(paletteLength) + 1

png(paste("results/heat_target_dn",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
dev.off()
}

if(trt=="ATC"){
  drugs_Bcell_sel<-names(Bcell_p_dn[[4]])[p.adjust(Bcell_p_dn[[4]], method="fdr")<=0.05]
  drugs_Tcell_sel<-names(Tcell_p_dn[[4]])[p.adjust(Tcell_p_dn[[4]], method="fdr")<=0.05]
  drugs_NKcell_sel<-names(NKcell_p_dn[[4]])[p.adjust(NKcell_p_dn[[4]], method="fdr")<=0.05]
  drugs_Myeloid_sel<-names(Myeloid_p_dn[[4]])[p.adjust(Myeloid_p_dn[[4]], method="fdr")<=0.05]
  drugs_allcells_dn<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel))
  
  drug_heat<-matrix(nrow=length(drugs_allcells_dn), ncol=4,
                    dimnames=list(c(drugs_allcells_dn), c("Bcell", "Tcell", "NKcell", "Myeloid")))
  
  drug_heat[,1]<-p.adjust(Bcell_p_dn[[4]], method="fdr")[drugs_allcells_dn]
  drug_heat[,2]<-p.adjust(Tcell_p_dn[[4]], method="fdr")[drugs_allcells_dn]
  drug_heat[,3]<-p.adjust(NKcell_p_dn[[4]], method="fdr")[drugs_allcells_dn]
  drug_heat[,4]<-p.adjust(Myeloid_p_dn[[4]], method="fdr")[drugs_allcells_dn]
  
  #drug_heat<-drug_heat[rowSums(drug_heat<0.05, na.rm=T)>2,]
  drug_heat[drug_heat<2.2*10^(-16)]<-2.2*10^(-16)
  drug_heat<-(-log10(drug_heat))
  
  if(nrow(drug_heat)>0){
    drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
    
    targets_dn<-rownames(drug_heat)
    
    
    
    paletteLength <- 50
    # use floor and ceiling to deal with even/odd length pallettelengths
    myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
    length(myBreaks) == length(paletteLength) + 1
    
    png(paste("results/heat_fourth_dn",trt,"enrich.png",sep="_"), res=600, 4000, 10000)
    print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
    dev.off()
  }
}
}

out<-list(drugs_up, moas_up, targets_up, drugs_dn, moas_dn, targets_dn)
return(out)
}


index_plot<-function(all_data, ids=drugs, class="drugs"){
  Bcell<-all_data[[1]]
  ind_up<-which(Bcell$fdr_q_nlog10>2 & Bcell$raw_cs>0)
  comp<-rep(NA, length(ind_up))
  if(class=="drugs"){
  comp[which(Bcell$pert_iname[ind_up] %in% ids)]<-Bcell$pert_iname[ind_up][which(Bcell$pert_iname[ind_up] %in% ids)]
  } else if (class=="moas"){
    comp[which(Bcell$moa[ind_up] %in% ids)]<-Bcell$moa[ind_up][which(Bcell$moa[ind_up] %in% ids)]
  }
  df<-data.frame(index=1:length(ind_up), drug=comp, cs=as.numeric(Bcell$norm_cs[ind_up]))
  df$drug<-factor(df$drug, levels=ids)
  plot<-ggplot(df, aes(x=index, y=cs))+geom_smooth()+
    geom_point(data=df[!is.na(df$drug),], aes(x=index, y=cs, colour=drug), size=2)+
    ylab("Connectivity score")+ggtitle("Bcells")+theme_bw()

  png(paste("results/Bcell",class, "index.png", sep="_"), res=600, 2500, 1500)
  print(plot)
  dev.off()
  
  Tcell<-all_data[[2]]
  ind_up<-which(Tcell$fdr_q_nlog10>2 & Tcell$raw_cs>0)
  comp<-rep(NA, length(ind_up))
  if(class=="drugs"){
    comp[which(Tcell$pert_iname[ind_up] %in% ids)]<-Tcell$pert_iname[ind_up][which(Tcell$pert_iname[ind_up] %in% ids)]
  } else if (class=="moas"){
    comp[which(Tcell$moa[ind_up] %in% ids)]<-Tcell$moa[ind_up][which(Tcell$moa[ind_up] %in% ids)]
  }
  df<-data.frame(index=1:length(ind_up), drug=comp, cs=as.numeric(Tcell$norm_cs[ind_up]))
  df$drug<-factor(df$drug, levels=ids)
  plot<-ggplot(df, aes(x=index, y=cs))+geom_smooth()+
    geom_point(data=df[!is.na(df$drug),], aes(x=index, y=cs, colour=drug), size=2)+
    ylab("Connectivity score")+ggtitle("Tcells")+theme_bw()
  png(paste("results/Tcell",class, "index.png", sep="_"), res=600, 2500, 1500)
  print(plot)
  dev.off()
  
  
  NKcell<-all_data[[3]]
  ind_up<-which(NKcell$fdr_q_nlog10>2 & NKcell$raw_cs>0)
  comp<-rep(NA, length(ind_up))
  if(class=="drugs"){
    comp[which(NKcell$pert_iname[ind_up] %in% ids)]<-NKcell$pert_iname[ind_up][which(NKcell$pert_iname[ind_up] %in% ids)]
  } else if (class=="moas"){
    comp[which(NKcell$moa[ind_up] %in% ids)]<-NKcell$moa[ind_up][which(NKcell$moa[ind_up] %in% ids)]
  }
  
  df<-data.frame(index=1:length(ind_up), drug=comp, cs=as.numeric(NKcell$norm_cs[ind_up]))
  df$drug<-factor(df$drug, levels=ids)
  plot<-ggplot(df, aes(x=index, y=cs))+geom_smooth()+
    geom_point(data=df[!is.na(df$drug),], aes(x=index, y=cs, colour=drug), size=2)+
    ylab("Connectivity score")+ggtitle("NKcells")+theme_bw()
  png(paste("results/NKcell",class, "index.png", sep="_"), res=600, 2500, 1500)
  print(plot)
  dev.off()
  
  Myeloid<-all_data[[4]]
  ind_up<-which(Myeloid$fdr_q_nlog10>2 & Myeloid$raw_cs>0)
  comp<-rep(NA, length(ind_up))
  if(class=="drugs"){
    comp[which(Myeloid$pert_iname[ind_up] %in% ids)]<-Myeloid$pert_iname[ind_up][which(Myeloid$pert_iname[ind_up] %in% ids)]
  } else if (class=="moas"){
    comp[which(Myeloid$moa[ind_up] %in% ids)]<-Myeloid$moa[ind_up][which(Myeloid$moa[ind_up] %in% ids)]
  }
  
  df<-data.frame(index=1:length(ind_up), drug=comp, cs=as.numeric(Myeloid$norm_cs[ind_up]))
  df$drug<-factor(df$drug, levels=ids)
  plot<-ggplot(df, aes(x=index, y=cs))+geom_smooth()+
    geom_point(data=df[!is.na(df$drug),], aes(x=index, y=cs, colour=drug), size=2)+
    ylab("Connectivity score")+ggtitle("Myeloid")+theme_bw()
  png(paste("results/Myeloid",class, "index.png", sep="_"), res=600, 2500, 1500)
  print(plot)
  dev.off()
  
}

library(ggplot2)



data_trt<-get_data(trt="trt_cp")
enrichment<-enrich_data(all_data=data_trt)
ids<-enrich_trt(all_enrich=enrichment, trt="trt_cp")

drugs<-ids[[1]]
index_plot(all_data=data_trt,ids=drugs, class="drugs")

moas<-ids[[2]]
index_plot(all_data=data_trt,ids=moas, class="moas")

###number of  hits
df<-data.frame(cell=rep(c("Bcells", "Tcells", "NKcells", "Myeloid"),3), number=c(length((data_trt[[1]]$pert_iname[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)])),
                                                                                               length((data_trt[[2]]$pert_iname[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)])),
                                                                                                             length((data_trt[[3]]$pert_iname[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)])),
                                                                                                                           length((data_trt[[4]]$pert_iname[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)])),
                                                                                 length((unlist(strsplit(data_trt[[1]]$moa[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[2]]$moa[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[3]]$moa[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[4]]$moa[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[1]]$target_name[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[2]]$target_name[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[3]]$target_name[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length((unlist(strsplit(data_trt[[4]]$target_name[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)], "\\|"))))),
               type=rep(c("drug", "moa", "target"),each=4)
)

png("results/Number_hits.png", res=600, 2500, 2500)
ggplot(df, aes(x=cell, fill=type, y=number))+geom_col(position="dodge", colour = "black")+theme_bw()+
  scale_fill_brewer(palette="Set2") 
dev.off()

df<-data.frame(cell=rep(c("Bcells", "Tcells", "NKcells", "Myeloid"),3), number=c(length(unique(data_trt[[1]]$pert_iname[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)])),
                                                                                 length(unique(data_trt[[2]]$pert_iname[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)])),
                                                                                 length(unique(data_trt[[3]]$pert_iname[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)])),
                                                                                 length(unique(data_trt[[4]]$pert_iname[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)])),
                                                                                 length(unique(unlist(strsplit(data_trt[[1]]$moa[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[2]]$moa[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[3]]$moa[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[4]]$moa[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[1]]$target_name[(data_trt[[1]]$raw_cs>0 & data_trt[[1]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[2]]$target_name[(data_trt[[2]]$raw_cs>0 & data_trt[[2]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[3]]$target_name[(data_trt[[3]]$raw_cs>0 & data_trt[[3]]$fdr_q_nlog10>2)], "\\|")))),
                                                                                 length(unique(unlist(strsplit(data_trt[[4]]$target_name[(data_trt[[4]]$raw_cs>0 & data_trt[[4]]$fdr_q_nlog10>2)], "\\|"))))),
               type=rep(c("drug", "moa", "target"),each=4)
)

png("results/Number_hits_unique.png", res=600, 2500, 2500)
ggplot(df, aes(x=cell, fill=type, y=number))+geom_col(position="dodge", colour = "black")+theme_bw()+
  scale_fill_brewer(palette="Set2") 
dev.off()

###number significant tests
df<-data.frame(cell=rep(c("Bcells", "Tcells", "NKcells", "Myeloid"),3), number=c(sum(p.adjust((enrichment[[1]][[1]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[2]][[1]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[3]][[1]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[4]][[1]]), "fdr")<0.05),
               sum(p.adjust((enrichment[[1]][[2]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[2]][[2]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[3]][[2]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[4]][[2]]), "fdr")<0.05),
               sum(p.adjust((enrichment[[1]][[3]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[2]][[3]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[3]][[3]]), "fdr")<0.05),
                       sum(p.adjust((enrichment[[4]][[3]]), "fdr")<0.05)),
               type=rep(c("drug", "moa", "target"),each=4)
)

png("results/Number_sig_tests.png", res=600, 2500, 2500)
ggplot(df, aes(x=cell, fill=type, y=number))+geom_col(position="dodge", colour = "black")+theme_bw()+
  scale_fill_brewer(palette="Set2") 
dev.off()

####Number of tests
Bcell<-read.csv("data/Bcell.gct", sep="\t", header=T)
Bcell<-Bcell[-1,]

df<-data.frame(type=c("drugs", "shRNA", "CRISPR", "OE"), 
               number=c(sum(Bcell$pert_type=="trt_cp"),
                        sum(Bcell$pert_type=="trt_sh.cgs"),
                        sum(Bcell$pert_type=="trt_xpr"),
                        sum(Bcell$pert_type=="trt_oe"))
)

png("results/Number_tests.png", res=600, 1500, 2500)
ggplot(df, aes(x=type, y=number))+geom_col(position="dodge", colour = "black")+theme_bw()+
  scale_fill_brewer(palette="Set2") 
dev.off()




