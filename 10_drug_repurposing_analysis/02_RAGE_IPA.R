setwd("working directory")


data<-read.csv("data/IPA_mild.gct", sep="\t", header=T)
data<-data[-1,]
data<-data[data$pert_type=="trt_cp",]

ind<-which(data$fdr_q_nlog10>2 & data$raw_cs>0)


drug_dn<-names(which(table(data$pert_iname[ind])>1))
drug_dn_p<-c()
for(drug in drug_dn){
  a<-sum(data$pert_iname[ind]==drug)
  c<-sum(data$pert_iname==drug)-a
  b<-length(ind)-sum(data$pert_iname[ind]==drug)
  d<-nrow(data)-length(unique(c(which(data$pert_iname==drug),ind )))
  drug_dn_p<-c(drug_dn_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
}
if(length(drug_dn)>0){
  names(drug_dn_p)<-drug_dn
}
names(drug_dn_p)[which(p.adjust(drug_dn_p, method="fdr")<0.05)]
mild<-p.adjust(drug_dn_p, method="fdr")[which(p.adjust(drug_dn_p, method="fdr")<0.05)]


rand<-function(data){
  
  ind<-sample(1:nrow(data), length(ind), replace = F)
  
  drug_Bcell_up<-names(which(table(data$pert_iname[ind])>1))
  drug_Bcell_up_p<-c()
  for(drug in drug_Bcell_up){
    a<-sum(data$pert_iname[ind]==drug)
    c<-sum(data$pert_iname==drug)-a
    b<-length(ind)-sum(data$pert_iname[ind]==drug)
    d<-nrow(data)-length(unique(c(which(data$pert_iname==drug),ind )))
    drug_Bcell_up_p<-c(drug_Bcell_up_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
  }
  if(length(drug_Bcell_up)>0){
    names(drug_Bcell_up_p)<-drug_Bcell_up
  }
  
  out<-p.adjust(drug_Bcell_up_p, method = "fdr")
  return(out)
}

set.seed(93748)
rand_drug<-list()
trials<-seq(1:10000)
parallel_enrich<-function(iter){
  print(iter)
  rand_drug<-rand(data=data)
  return(rand_drug)
}

library(parallel)
RAGE_drug_rand <- mclapply(trials, parallel_enrich, mc.cores = 1)
save(RAGE_drug_rand, file="RAGE_drug_rand.RData")

random_drugs<-table(names(unlist(RAGE_drug_rand))[unlist(RAGE_drug_rand)<0.05])
#random_drugs[names(mild)]/10000

IPA<-read.xlsx("data/PredictionIpa.xlsx",1)
targets<-unlist(strsplit(data$target_name, "\\|"))
targets_sel<-unlist(strsplit(data$target_name[ind], "\\|"))

target_int<-intersect(targets, IPA$Symbol)
drugs_sel<-c()
for(i in 1:length(target_int)){
drugs_sel<-c(drugs_sel, unique(data$pert_iname[grep(target_int[i], data$target_name)]))
}

intersect(drugs_sel, names(mild))

data<-read.csv("data/IPA_severe.gct", sep="\t", header=T)
data<-data[-1,]
data<-data[data$pert_type=="trt_cp",]

ind<-which(data$fdr_q_nlog10>2 & data$raw_cs>0)


drug_dn<-names(which(table(data$pert_iname[ind])>1))
drug_dn_p<-c()
for(drug in drug_dn){
  a<-sum(data$pert_iname[ind]==drug)
  c<-sum(data$pert_iname==drug)-a
  b<-length(ind)-sum(data$pert_iname[ind]==drug)
  d<-nrow(data)-length(unique(c(which(data$pert_iname==drug),ind )))
  drug_dn_p<-c(drug_dn_p, fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")[[1]])
}
if(length(drug_dn)>0){
  names(drug_dn_p)<-drug_dn
}
names(drug_dn_p)[which(p.adjust(drug_dn_p, method="fdr")<0.05)]
severe<-p.adjust(drug_dn_p, method="fdr")[which(p.adjust(drug_dn_p, method="fdr")<0.05)]



inboth<-intersect(names(mild), names(severe))
cbind(mild[inboth], severe[inboth])
write.csv(mild, "mild.csv")
write.csv(severe, "severe.csv")

