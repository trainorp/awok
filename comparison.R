############ Prereqs ############
# ptm<-proc.time()
options(stringsAsFactors=FALSE)
library(randomForest)
library(e1071)
library(nnet)
library(glmnet)
library(cvTools)
library(dplyr)
library(tidyr)
library(doParallel)
library(ggplot2)
library(gridExtra)

oldPar<-par()

setwd("~/gdrive/AthroMetab/WoAC")

## Do not run:
setwd("~/gdrive/AthroMetab/Data")
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

setwd("~/gdrive/AthroMetab/WoAC")

########### Random Forest importance ###########
# cl<-makeCluster(10)
# registerDoParallel(cl)

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

# print(proc.time()-ptm)
# stopCluster(cl)

save.image(file="compare.RData")

## End do not run:
load(file="compare.RData")

# RF plot:
plot(bfe$iter,bfe$err,type="l")

p1<-ggplot(bfe %>% filter(is.finite(err)),aes(x=iter,y=err))+geom_line()+theme_bw()+
  xlab("Features (metabolites) eliminated")+ylab("Cross-entropy error (Repeated 10-fold CV)")

########### Lasso / Elastic net ###########
# List of elastic nets:
tGLMlist<-list()
for(alpha in seq(.1,1,.1)){
  tGLM<-glmnet(x=as.matrix(metab[,names(metab)!="group"]),alpha=alpha,
               y=metab[,"group"],family="multinomial",lambda=10**(-seq(0.44,6,.02)))
  tGLMlist[[as.character(alpha)]]<-tGLM
}
# CV version:
tGLMcvlist<-list()
for(alpha in seq(.1,1,.1)){
  tGLMcv<-cv.glmnet(x=as.matrix(metab[,names(metab)!="group"]),alpha=alpha,
               y=metab[,"group"],family="multinomial",lambda=10**(-seq(0.44,6,.02)),
               type.measure="deviance")
  tGLMcvlist[[as.character(alpha)]]<-tGLMcv
}

tGLM<-tGLMlist$`0.9`
tGLMdf<-data.frame(lambda=tGLM$lambda,metabs=NA,l1Norm=NA)
for(i in 1:nrow(tGLMdf)){
  tGLMcoef<-coef(tGLM,s=tGLMdf$lambda[i])
  tGLMcoef<-cbind(as.matrix(tGLMcoef$sCAD),as.matrix(tGLMcoef$Type.1.MI),
                  as.matrix(tGLMcoef$Type.2.MI))
  colnames(tGLMcoef)<-c("sCAD","Type.1.MI","Type.2.MI")
  if(!is.null(nrow(tGLMcoef[apply(tGLMcoef,1,FUN=function(x) sum(as.integer(x==0))<3),]))){
    tGLMdf$metabs[i]<-nrow(tGLMcoef[apply(tGLMcoef,1,FUN=function(x) sum(as.integer(x==0))<3),])-1
    tGLMdf$l1Norm[i]<-sum(abs(tGLMcoef))
  }
}

# Plot of L1 norm vs features selected
tGLMdfLong<-tGLMdf %>% gather(key="Measure",value="value",-lambda)
p2<-ggplot(tGLMdfLong %>% filter(!is.na(value)),aes(x=lambda,y=value,color=Measure))+geom_line()+theme_bw()+
  xlab(expression(lambda))+ylab("Value")+
  scale_color_manual(values=c("darkred","navyblue"),
    labels=c("L1 Norm","Features (metabolites)\nselected"))

GLMerr<-data.frame(lambda=tGLMcv$lambda,err=log2(exp(tGLMcv$cvm))/2)
tGLMdf<-tGLMdf %>% full_join(GLMerr)
p3<-ggplot(tGLMdf ,aes(x=metabs,y=err))+geom_line()+theme_bw()+
  xlab("Features (metabolites) selected")+ylab("Cross-entropy error (Repeated 10-fold CV)")
p4<-ggplot(tGLMdf,aes(x=lambda,y=err))+geom_line()+theme_bw()+
  xlab(expression(lambda))+ylab("Cross-entropy error (Repeated 10-fold CV)")

grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
# Check
tGLMcv<-tGLMcvlist$`0.9`
set.seed(3)
cvF1<-cvFolds(n=nrow(metab),K=10,R=10)
elasticp<-foreach(j=1:cvF1$R, .combine="rbind",.packages="glmnet") %dopar% {
  elasticp<-data.frame()
  # K-folds:
  for(k in 1:cvF1$K){
    cvF1Df<-data.frame(ind=cvF1$subsets[,j],fold=cvF1$which)
    elastic1<-glmnet(x=as.matrix(metab[cvF1Df$ind[cvF1Df$fold!=k],names(metab)!="group"]),
                y=metab[cvF1Df$ind[cvF1Df$fold!=k],"group"],alpha=.9,
                family="multinomial",lambda=.2)
    elastic1p<-predict(elastic1,newx=as.matrix(metab[cvF1Df$ind[cvF1Df$fold==k],names(metab)!="group"]),
                  type="response")
    elastic1p<-as.data.frame(elastic1p)
    elastic1p$ptid<-rownames(elastic1p)
    elasticp<-rbind(elasticp,elastic1p)
  }
  elasticp
}
elasticp$group<-metab$group[match(elasticp$ptid,rownames(metab))]
elasticp$ptid<-NULL
colnames(elasticp)<-gsub(".s0","",colnames(elasticp))
elasticp$loss<--log(
  as.numeric(as.matrix(elasticp)[1:nrow(elasticp) + nrow(elasticp) * (match(elasticp$group,colnames(elasticp)) - 1)]))
mean(elasticp$loss)
