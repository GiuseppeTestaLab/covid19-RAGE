library(pheatmap)
setwd("~/Work/Projects/Covidiamo/cMAP/covidiamo_drugs")
source("covidiamo_functions_RAGE.R")

data_trt<-get_data(trt="trt_cp")
gene_ind<-get_gene_ind(all_data=data_trt)
enrichment<-enrich_data(all_data=data_trt, ind=gene_ind)
ids<-enrich_trt(all_enrich=enrichment, thr_up=1,  trt="trt_cp", to_plot=T)


Myeloid_p_up<-enrichment[[4]]
drugs_Myeloid_sel<-names(Myeloid_p_up[[1]])[p.adjust(Myeloid_p_up[[1]], method="fdr")<=0.05]


Bcell_p_up<-enrichment[[1]]
Tcell_p_up<-enrichment[[2]]
NKcell_p_up<-enrichment[[3]]
Myeloid_p_up<-enrichment[[4]]
CP_p_up<-enrichment[[5]]

drug_heat<-matrix(nrow=length(drugs_Myeloid_sel), ncol=5,
                  dimnames=list(c(drugs_Myeloid_sel), c("Bcell", "Tcell", "NKcell", "Myeloid", "CP")))

drug_heat[,1]<-p.adjust(Bcell_p_up[[1]], method="fdr")[drugs_Myeloid_sel]
drug_heat[,2]<-p.adjust(Tcell_p_up[[1]], method="fdr")[drugs_Myeloid_sel]
drug_heat[,3]<-p.adjust(NKcell_p_up[[1]], method="fdr")[drugs_Myeloid_sel]
drug_heat[,4]<-p.adjust(Myeloid_p_up[[1]], method="fdr")[drugs_Myeloid_sel]
drug_heat[,5]<-p.adjust(CP_p_up[[1]], method="fdr")[drugs_Myeloid_sel]
drug_heat<-(-log10(drug_heat))
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
    
    png(paste("results/heat_drugs_up_cp_enrich_RAGE_myel.png",sep="_"), res=600, 2000, 16000)
    print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
    dev.off()
 