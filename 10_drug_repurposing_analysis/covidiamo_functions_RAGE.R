##COVIDIAMO randomizations
#reshuffle labels 10000
#how many drugs/moas/targets pass the selection filters?
#how many times the specifically selected drugs/moas/targets pass the selection filters?
get_data<-function(trt){
  Bcell<-read.csv("data/B_RAGE.gct", sep="\t", header=T)
  Bcell<-Bcell[-1,]
  Bcell<-Bcell[Bcell$pert_type==trt,]
  
  
  Tcell<-read.csv("data/T_RAGE.gct",  sep="\t", header=T)
  Tcell<-Tcell[-1,]
  Tcell<-Tcell[Tcell$pert_type==trt,]
  
  
  NKcell<-read.csv("data/NK_RAGE.gct", sep="\t", header=T)
  NKcell<-NKcell[-1,]
  NKcell<-NKcell[NKcell$pert_type==trt,]
  
  
  Myeloid<-read.csv("data/Myeloid_RAGE.gct", sep="\t", header=T)
  Myeloid<-Myeloid[-1,]
  Myeloid<-Myeloid[Myeloid$pert_type==trt,]
  
  CP<-read.csv("data/CP_RAGE.gct", sep="\t", header=T)
  CP<-CP[-1,]
  CP<-CP[CP$pert_type==trt,]
  
  all_data<-list(Bcell, Tcell, NKcell, Myeloid, CP)
  return(all_data)
}


get_gene_ind<-function(all_data){
  ind<-list()
  Bcell<-all_data[[1]]
  ind[[1]]<-which(Bcell$fdr_q_nlog10>2 & Bcell$raw_cs>0)
  ind[[2]]<-which(Bcell$fdr_q_nlog10>2 & Bcell$raw_cs<0)
  Tcell<-all_data[[2]]
  ind[[3]]<-which(Tcell$fdr_q_nlog10>2 & Tcell$raw_cs>0)
  ind[[4]]<-which(Tcell$fdr_q_nlog10>2 & Tcell$raw_cs<0)
  NKcell<-all_data[[3]]
  ind[[5]]<-which(NKcell$fdr_q_nlog10>2 & NKcell$raw_cs>0)
  ind[[6]]<-which(NKcell$fdr_q_nlog10>2 & NKcell$raw_cs<0)
  Myeloid<-all_data[[4]]
  ind[[7]]<-which(Myeloid$fdr_q_nlog10>2 & Myeloid$raw_cs>0)
  ind[[8]]<-which(Myeloid$fdr_q_nlog10>2 & Myeloid$raw_cs<0)
  CP<-all_data[[5]]
  ind[[9]]<-which(CP$fdr_q_nlog10>2 & CP$raw_cs>0)
  ind[[10]]<-which(CP$fdr_q_nlog10>2 & CP$raw_cs<0)
  names(ind)<-c("Bcell_up", "Bcell_dn",
                "Tcell_up", "Tcell_dn",
                "NKcell_up", "NKcell_dn",
                "Myeloid_up", "Myeloid_dn",
                "CP_up", "CP_dn")
  return(ind)
  
}


get_gene_ind_rand<-function(all_data, ind){
  ind_rand<-list()
  Bcell<-all_data[[1]]
  ind_rand[[1]]<-sample(1:nrow(Bcell), length(ind[[1]]), replace = F)
  ind_rand[[2]]<-sample(1:nrow(Bcell), length(ind[[2]]), replace = F)
  Tcell<-all_data[[2]]
  ind_rand[[3]]<-sample(1:nrow(Tcell), length(ind[[3]]), replace = F)
  ind_rand[[4]]<-sample(1:nrow(Tcell), length(ind[[4]]), replace = F)
  NKcell<-all_data[[3]]
  ind_rand[[5]]<-sample(1:nrow(NKcell), length(ind[[5]]), replace = F)
  ind_rand[[6]]<-sample(1:nrow(NKcell), length(ind[[6]]), replace = F)
  Myeloid<-all_data[[4]]
  ind_rand[[7]]<-sample(1:nrow(Myeloid), length(ind[[7]]), replace = F)
  ind_rand[[8]]<-sample(1:nrow(Myeloid), length(ind[[8]]), replace = F)
  CP<-all_data[[5]]
  ind_rand[[9]]<-sample(1:nrow(CP), length(ind[[9]]), replace = F)
  ind_rand[[10]]<-sample(1:nrow(CP), length(ind[[10]]), replace = F)
  
  names(ind_rand)<-c("Bcell_up", "Bcell_dn",
                "Tcell_up", "Tcell_dn",
                "NKcell_up", "NKcell_dn",
                "Myeloid_up", "Myeloid_dn")
  return(ind_rand)
  
}

enrich<-function(cMAP, ind){
  drug_Bcell_up<-names(which(table(cMAP$pert_iname[ind])>1))
  drug_Bcell_up_p<-c()
  for(drug in drug_Bcell_up){
    a<-sum(cMAP$pert_iname[ind]==drug)
    c<-sum(cMAP$pert_iname==drug)-a
    b<-length(ind)-sum(cMAP$pert_iname[ind]==drug)
    d<-nrow(cMAP)-length(unique(c(which(cMAP$pert_iname==drug),ind )))
    drug_Bcell_up_p<-c(drug_Bcell_up_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
  }
  if(length(drug_Bcell_up)>0){
  names(drug_Bcell_up_p)<-drug_Bcell_up
  }

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
  if(length(drug_Bcell_up)>0){
  names(moa_Bcell_up_p)<-moa_Bcell_up
  }
  
  
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
  if(length(target_Bcell_up)>0){
  names(target_Bcell_up_p)<-target_Bcell_up
  }
  pvals<-list(drug=drug_Bcell_up_p, moa=moa_Bcell_up_p, target=target_Bcell_up_p)
  return(pvals)
}
  
enrich_data<-function(all_data, ind){
   
  Bcell<-all_data[[1]]
  Bcell_p_up<-enrich(cMAP=Bcell, ind=ind[[1]])
  Bcell_p_dn<-enrich(cMAP=Bcell, ind=ind[[2]])
 
  Tcell<-all_data[[2]]
  Tcell_p_up<-enrich(cMAP=Tcell, ind=ind[[3]])
  Tcell_p_dn<-enrich(cMAP=Tcell, ind=ind[[4]])
  
  NKcell<-all_data[[3]]
  NKcell_p_up<-enrich(cMAP=NKcell, ind=ind[[5]])
  NKcell_p_dn<-enrich(cMAP=NKcell, ind=ind[[6]])
  
  
  Myeloid<-all_data[[4]]
  Myeloid_p_up<-enrich(cMAP=Myeloid, ind=ind[[7]])
  Myeloid_p_dn<-enrich(cMAP=Myeloid, ind=ind[[8]])
  
  CP<-all_data[[5]]
  CP_p_up<-enrich(cMAP=CP, ind=ind[[9]])
  CP_p_dn<-enrich(cMAP=CP, ind=ind[[10]])
  
  all_enrich<-list(Bcell_p_up, Tcell_p_up, NKcell_p_up, Myeloid_p_up,CP_p_up, 
                   Bcell_p_dn, Tcell_p_dn, NKcell_p_dn, Myeloid_p_dn,
                   CP_p_dn)
  return(all_enrich)
}

enrich_trt<-function(all_enrich, trt, to_plot=F, thr_up=2, thr_dn=1){
  Bcell_p_up<-all_enrich[[1]]
  Tcell_p_up<-all_enrich[[2]]
  NKcell_p_up<-all_enrich[[3]]
  Myeloid_p_up<-all_enrich[[4]]
  CP_p_up<-all_enrich[[5]]
  Bcell_p_dn<-all_enrich[[6]]
  Tcell_p_dn<-all_enrich[[7]]
  NKcell_p_dn<-all_enrich[[8]]
  Myeloid_p_dn<-all_enrich[[9]]
  CP_p_dn<-all_enrich[[10]]
  
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
  drugs_CP_sel<-names(CP_p_up[[1]])[p.adjust(CP_p_up[[1]], method="fdr")<=0.05]
  drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel, drugs_CP_sel))
  
  drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=5,
                    dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid", "CP")))
  
  drug_heat[,1]<-p.adjust(Bcell_p_up[[1]], method="fdr")[drugs_allcells_up]
  drug_heat[,2]<-p.adjust(Tcell_p_up[[1]], method="fdr")[drugs_allcells_up]
  drug_heat[,3]<-p.adjust(NKcell_p_up[[1]], method="fdr")[drugs_allcells_up]
  drug_heat[,4]<-p.adjust(Myeloid_p_up[[1]], method="fdr")[drugs_allcells_up]
  drug_heat[,5]<-p.adjust(CP_p_up[[1]], method="fdr")[drugs_allcells_up]
  
  ind_heat<-which(rowSums(drug_heat<0.05, na.rm=T)>thr_up)
  if(!is.null(ind_heat)&length(ind_heat)==1){
    drugs_up<-names(ind_heat)
  } else if (!is.null(ind_heat)&length(ind_heat)>1){
    
  drug_heat<-drug_heat[ind_heat,]
  drug_heat[which(drug_heat<2.2*10^(-16))]<-2.2*10^(-16)
  drug_heat<-(-log10(drug_heat))
  
    drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
    drugs_up<-rownames(drug_heat)
    if(to_plot){
    paletteLength <- 50
    # use floor and ceiling to deal with even/odd length pallettelengths
    myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(unlist(drug_heat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=floor(paletteLength/2)))
    length(myBreaks) == length(paletteLength) + 1
    
    png(paste("results_RAGE/heat_drugs_up",trt,"enrich_RAGE.png",sep="_"), res=600, 2000, 8000)
    print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
    dev.off()
  }
  }
  

    
    #################
    ##moa
    ###################
    
    drugs_Bcell_sel<-names(Bcell_p_up[[2]])[p.adjust(Bcell_p_up[[2]], method="fdr")<=0.05]
    drugs_Tcell_sel<-names(Tcell_p_up[[2]])[p.adjust(Tcell_p_up[[2]], method="fdr")<=0.05]
    drugs_NKcell_sel<-names(NKcell_p_up[[2]])[p.adjust(NKcell_p_up[[2]], method="fdr")<=0.05]
    drugs_Myeloid_sel<-names(Myeloid_p_up[[2]])[p.adjust(Myeloid_p_up[[2]], method="fdr")<=0.05]
    drugs_CP_sel<-names(CP_p_up[[2]])[p.adjust(CP_p_up[[2]], method="fdr")<=0.05]
    drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel, drugs_CP_sel))
    
    drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=5,
                      dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid", "CP")))
    
    drug_heat[,1]<-p.adjust(Bcell_p_up[[2]], method="fdr")[drugs_allcells_up]
    drug_heat[,2]<-p.adjust(Tcell_p_up[[2]], method="fdr")[drugs_allcells_up]
    drug_heat[,3]<-p.adjust(NKcell_p_up[[2]], method="fdr")[drugs_allcells_up]
    drug_heat[,4]<-p.adjust(Myeloid_p_up[[2]], method="fdr")[drugs_allcells_up]
    drug_heat[,5]<-p.adjust(CP_p_up[[2]], method="fdr")[drugs_allcells_up]
    
    ind_heat<-which(rowSums(drug_heat<0.05, na.rm=T)>1)
    if(!is.null(ind_heat)&length(ind_heat)==1){
      moas_up<-names(ind_heat)
    } else if (!is.null(ind_heat)&length(ind_heat)>1){
      
    drug_heat<-drug_heat[ind_heat,]
    drug_heat[which(drug_heat<2.2*10^(-16))]<-2.2*10^(-16)
    drug_heat<-(-log10(drug_heat))
    
      drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
      moas_up<-rownames(drug_heat)
      if(to_plot){
      
      paletteLength <- 50
      # use floor and ceiling to deal with even/odd length pallettelengths
      myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
      # length(breaks) == length(paletteLength) + 1
      # use floor and ceiling to deal with even/odd length pallettelengths
      myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
      length(myBreaks) == length(paletteLength) + 1
      
      png(paste("results_RAGE/heat_moa_up",trt,"enrich_RAGE.png",sep="_"), res=600, 4000, 10000)
      print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
      dev.off()
      }
    }
    
    #################
    ##targets
    ###################
    
    drugs_Bcell_sel<-names(Bcell_p_up[[3]])[p.adjust(Bcell_p_up[[3]], method="fdr")<=0.05]
    drugs_Tcell_sel<-names(Tcell_p_up[[3]])[p.adjust(Tcell_p_up[[3]], method="fdr")<=0.05]
    drugs_NKcell_sel<-names(NKcell_p_up[[3]])[p.adjust(NKcell_p_up[[3]], method="fdr")<=0.05]
    drugs_Myeloid_sel<-names(Myeloid_p_up[[3]])[p.adjust(Myeloid_p_up[[3]], method="fdr")<=0.05]
    drugs_CP_sel<-names(CP_p_up[[3]])[p.adjust(CP_p_up[[3]], method="fdr")<=0.05]
    drugs_allcells_up<-unique(c(drugs_Bcell_sel, drugs_Tcell_sel, drugs_NKcell_sel, drugs_Myeloid_sel, drugs_CP_sel))
    
    drug_heat<-matrix(nrow=length(drugs_allcells_up), ncol=5,
                      dimnames=list(c(drugs_allcells_up), c("Bcell", "Tcell", "NKcell", "Myeloid", "CP")))
    
    drug_heat[,1]<-p.adjust(Bcell_p_up[[3]], method="fdr")[drugs_allcells_up]
    drug_heat[,2]<-p.adjust(Tcell_p_up[[3]], method="fdr")[drugs_allcells_up]
    drug_heat[,3]<-p.adjust(NKcell_p_up[[3]], method="fdr")[drugs_allcells_up]
    drug_heat[,4]<-p.adjust(Myeloid_p_up[[3]], method="fdr")[drugs_allcells_up]
    drug_heat[,5]<-p.adjust(CP_p_up[[3]], method="fdr")[drugs_allcells_up]
    
    ind_heat<-which(rowSums(drug_heat<0.05, na.rm=T)>1)
    if(!is.null(ind_heat)&length(ind_heat)==1){
      targets_up<-names(ind_heat)
    } else if (!is.null(ind_heat)&length(ind_heat)>1){
    drug_heat<-drug_heat[ind_heat,]
    drug_heat[which(drug_heat<2.2*10^(-16))]<-2.2*10^(-16)
    drug_heat<-(-log10(drug_heat))
    
    
      
      drug_heat<-drug_heat[order(rowMeans(drug_heat), decreasing=T),]
      
      targets_up<-rownames(drug_heat)
      if(to_plot){
      
      paletteLength <- 50
      # use floor and ceiling to deal with even/odd length pallettelengths
      myColor <- colorRampPalette(c( "white", "#D73027"))(paletteLength)
      # length(breaks) == length(paletteLength) + 1
      # use floor and ceiling to deal with even/odd length pallettelengths
      myBreaks <- c(0, seq(max(unlist(drug_heat), na.rm=T)/paletteLength, max(unlist(drug_heat), na.rm=T), length.out=paletteLength-1))
      length(myBreaks) == length(paletteLength) + 1
      
      png(paste("results_RAGE/heat_target_up",trt,"enrich_RAGE.png",sep="_"), res=600, 4000, 10000)
      print(pheatmap(drug_heat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F))
      dev.off()
      }
    } 
    
    
  out<-list(drugs_up, moas_up, targets_up, drugs_dn, moas_dn, targets_dn)
  return(out)
}



