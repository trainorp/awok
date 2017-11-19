############ Prereqs ############
ptm<-proc.time()
options(stringsAsFactors=FALSE)
library(randomForest)
library(e1071)
library(nnet)
library(cvTools)
library(dplyr)
library(doParallel)

# oldPar<-par()

# setwd("~/gdrive/AthroMetab/Data")
spec<-read.csv("AthroACSRawSpectra.csv")
key<-read.csv("metabolite_key2.csv")
key[grepl("X - ",key$biochemical),]$biochemical<-key[grepl("X - ",key$biochemical),]$Unknown.Name

key<-spec %>% left_join(key,c("comp_id"="MEBID"))

getMZ<-function(x)
{
  x<-strsplit(x," ")[[1]]
  return(unlist(strsplit(x[grepl(":100",x)],":"))[1])
}
key$peakMZ<-unlist(lapply(lapply(key$spectra,getMZ),function(x) ifelse(is.null(x),NA,x)))

options(stringsAsFactors=FALSE)

############ Import ############
metab<-read.csv("scaled.csv")
pheno<-metab %>% filter(timepoint=="T0") %>% select(ptid,group)
metab<-metab %>% filter(timepoint=="T0") %>% select(-ptid,-group,-timepoint) 

#Entropy filter & log-transform
metab<-metab[,apply(metab,2,function(x) table(x)[1]<26)] %>% log2()
rownames(metab)<-pheno$ptid
pheno$group2<-ifelse(pheno$group=="Type 1 MI","Thromb. MI",
                     ifelse(pheno$group=="Type 2 MI","Non-Thromb. MI","sCAD"))
pheno$ptid2<-NA
set.seed(1)
pheno$ptid2[pheno$group=="sCAD"]<-sample(1:length(pheno$ptid2[pheno$group=="sCAD"]))
set.seed(2)
pheno$ptid2[pheno$group=="Type 1 MI"]<-sample(1:length(pheno$ptid2[pheno$group=="Type 1 MI"]))
set.seed(3)
pheno$ptid2[pheno$group=="Type 2 MI"]<-sample(1:length(pheno$ptid2[pheno$group=="Type 2 MI"]))

# Add phenotype back:
metab<-cbind(group=pheno$group,metab)
metab$group<-factor(make.names(metab$group))

# setwd("~/gdrive/AthroMetab/WoAC")

########### Random Forest importance ###########
cl<-makeCluster(10)
registerDoParallel(cl)

metab1<-metab
bfe<-data.frame(iter=1:(ncol(metab1)-1))
bfe$err<-NA
bfe$elim<-""

for(i in 1:nrow(bfe)){
  set.seed(i)
  cvF1<-cvFolds(n=nrow(metab1),K=10,R=10)
  
  # Repeats for repeated K-fold CV:
  rfp<-foreach(j=1:cvF1$R, .combine="rbind",.packages="randomForest") %dopar% {
    rfp<-data.frame()
    # K-folds:
    for(k in 1:cvF1$K){
      cvF1Df<-data.frame(ind=cvF1$subsets[,j],fold=cvF1$which)
      rf1<-randomForest(group~.,data=metab1[cvF1Df$ind[cvF1Df$fold!=k],],
                        ntree=500,importance=TRUE)
      rf1p<-predict(rf1,newdata=metab1[cvF1Df$ind[cvF1Df$fold==k],],type="prob")
      rf1p<-as.data.frame(rf1p)
      rf1p$ptid<-rownames(rf1p)
      rfp<-rbind(rfp,rf1p)
    }
    rfp
  }
  
  rfp$group<-metab1$group[match(rfp$ptid,rownames(metab1))]
  rfp$ptid<-NULL
  rfp$loss<--log2(
    as.numeric(as.matrix(rfp)[1:nrow(rfp) + nrow(rfp) * (match(rfp$group,colnames(rfp)) - 1)]))
  bfe$err[i]<-mean(rfp$loss)
  
  rf1<-randomForest(group~.,data=metab1,ntree=1000,importance=TRUE)
  imp1<-as.data.frame(importance(rf1))
  imp1$metab1<-rownames(imp1)
  metab1<-metab1[,names(metab1)!=imp1$metab1[which.min(imp1$MeanDecreaseAccuracy)]]
  bfe$elim[i]<-imp1$metab1[which.min(imp1$MeanDecreaseAccuracy)]
  print(i)
}

print(proc.time()-ptm)
stopCluster(cl)

save.image(file="compare.RData")
