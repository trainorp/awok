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

########### Random Forest importance ###########
# cl<-makeCluster(10)
# registerDoParallel(cl)

rF0<-randomForest(group~.,data=metab,ntree=10000,importance=TRUE)

## BFE:
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

# RF plot:
plot(bfe$iter,bfe$err,type="l")

p1<-ggplot(bfe %>% filter(is.finite(err)),aes(x=iter,y=err))+geom_line()+theme_bw()+
  xlab("Features (metabolites) eliminated")+ylab("Cross-entropy error (Repeated 10-fold CV)")

########### Lasso / Elastic net ###########
# List of elastic nets:
tGLMlist<-list()
for(alpha in seq(.5,1,.1)){
  tGLM<-glmnet(x=as.matrix(metab[,names(metab)!="group"]),alpha=alpha,
               y=metab[,"group"],family="multinomial",lambda=10**(-seq(0.2,6,.001)))
  
  tGLMdf<-data.frame(lambda=tGLM$lambda,metabs=NA,l1Norm=NA)
  varList<-list()
  for(i in 1:nrow(tGLMdf)){
    tGLMcoef<-coef(tGLM,s=tGLMdf$lambda[i])
    tGLMcoef<-cbind(as.matrix(tGLMcoef$sCAD),as.matrix(tGLMcoef$Type.1.MI),
                    as.matrix(tGLMcoef$Type.2.MI))
    colnames(tGLMcoef)<-c("sCAD","Type.1.MI","Type.2.MI")
    if(!is.null(nrow(tGLMcoef[apply(tGLMcoef,1,FUN=function(x) sum(as.integer(x==0))<3),]))){
      varList[[i]]<-paste(rownames(tGLMcoef[apply(tGLMcoef,1,FUN=function(x) sum(as.integer(x==0))<3),]))
      tGLMdf$metabs[i]<-nrow(tGLMcoef[apply(tGLMcoef,1,FUN=function(x) sum(as.integer(x==0))<3),])-1
      tGLMdf$l1Norm[i]<-sum(abs(tGLMcoef))
    }
  }
  tGLM$nMetabDf<-tGLMdf
  tGLM$varList<-varList
  
  tGLMlist[[as.character(alpha)]]<-tGLM
}
# CV version:
tGLMcvlist<-list()
for(alpha in seq(.5,1,.1)){
  tGLMcv<-cv.glmnet(x=as.matrix(metab[,names(metab)!="group"]),alpha=alpha,
               y=metab[,"group"],family="multinomial",lambda=10**(-seq(0.2,6,.001)),
               type.measure="deviance",nfolds=10)
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
    labels=c(expression(paste(L[1]," Norm")),"Features (metabolites) selected"))+
  ggtitle(expression(paste("(C) Norm and selection over ",lambda," path")))+
  theme(legend.text.align = 0)

# Plot of error vs features
GLMerr<-data.frame(lambda=tGLMcv$lambda,err=log2(exp(tGLMcv$cvm))/2)
tGLMdf<-tGLMdf %>% full_join(GLMerr)
p3<-ggplot(tGLMdf ,aes(x=metabs,y=err))+geom_line()+theme_bw()+
  xlab("Features (metabolites) selected")+ylab("Cross-entropy error")+
  ggtitle("(B) Error versus # features")

# Plot of error vs lambda
p4<-ggplot(tGLMdf,aes(x=lambda,y=err))+geom_line()+theme_bw()+
  xlab(expression(lambda))+ylab("Cross-entropy error")+xlim(0,.35)+
  ggtitle(expression(paste("(A) Error over ",lambda," path")))

lm1<-matrix(c(3,3,3,3,3,3,3,2,2,2,2,2,2,2,NA,1,1,1,1,1,1,1,1,1,1,1,1,NA),nrow=2,byrow=TRUE)
# png(file="elasticNet.png",height=5,width=6,units="in",res=600)
grid.arrange(p2,p3,p4,layout_matrix=lm1)
# dev.off()

########### WoAC ###########
load("outs.RData")
weights<-c()
for(i in 1:length(outs)){
  nams<-outs[[i]]$vars
  experts<-outs[[i]]$pop@varInclude[order(outs[[i]]$pop@cost)[1]]
  experts<-sapply(experts,as.numeric)
  rownames(experts)<-nams
  blanks<-matrix(0,nrow=ncol(metab)-1,ncol=1)
  rownames(blanks)<-colnames(metab)[colnames(metab)!="group"]
  blanks[match(rownames(experts),rownames(blanks)),]<-experts
  weights<-cbind(weights,apply(blanks,1,mean))
}
key<-key %>% arrange(as.numeric(gsub("M","",id)))
weights<-apply(weights,1,function(x) mean(x))
key$weights<-0
key$weights[match(names(weights),key$id)]<-weights

save.image(file="compare.RData")

## End do not run:
setwd("~/gdrive/AthroMetab/WoAC")
load(file="compare.RData")

############ Selection ############
selDf<-key[,names(key)%in%c("id","biochemical","weights")]

# RF importance:
rF0Imp<-as.data.frame(importance(rF0,type="1"))
rF0Imp$id<-rownames(rF0Imp)
rF0Imp<-rF0Imp %>% arrange(desc(MeanDecreaseAccuracy))
rF0Imp$rFImp<-1:nrow(rF0Imp)
rF0Imp$MeanDecreaseAccuracy<-NULL

# RF-BFE importance:
rFBFE<-bfe
rFBFE<-rFBFE %>% arrange(desc(iter))
rFBFE$rFBFEImp<-1:max(rFBFE$iter)
rFBFE<-rFBFE %>% select(elim,rFBFEImp)
names(rFBFE)[names(rFBFE)=="elim"]<-"id"

# WoAC importance:
WoACImp<-key %>% select(id,weights)
WoACImp<-WoACImp %>% arrange(desc(weights))
WoACImp$WoACImp<-1:nrow(WoACImp)
WoACImp$weights<-NULL

selDf<-selDf %>% left_join(rF0Imp) %>% left_join(rFBFE) %>% left_join(WoACImp)
lasso<-tGLMlist[["1"]]

# 3 Metabolite:
selDf$lasso3<-0
lamb<-min(lasso$nMetabDf$lambda[!is.na(lasso$nMetabDf$metabs) & lasso$nMetabDf$metabs==3])
selDf$lasso3[selDf$id %in% lasso$varList[which(lasso$nMetabDf$lambda==lamb)][[1]]]<-1

# 5 Metabolite:
selDf$lasso5<-0
lamb<-min(lasso$nMetabDf$lambda[!is.na(lasso$nMetabDf$metabs) & lasso$nMetabDf$metabs==5])
selDf$lasso5[selDf$id %in% lasso$varList[which(lasso$nMetabDf$lambda==lamb)][[1]]]<-1

# 7 Metabolite:
selDf$lasso7<-0
lamb<-min(lasso$nMetabDf$lambda[!is.na(lasso$nMetabDf$metabs) & lasso$nMetabDf$metabs==7])
selDf$lasso7[selDf$id %in% lasso$varList[which(lasso$nMetabDf$lambda==lamb)][[1]]]<-1

# 10 Metabolite:
selDf$lasso10<-0
lamb<-min(lasso$nMetabDf$lambda[!is.na(lasso$nMetabDf$metabs) & lasso$nMetabDf$metabs==10])
selDf$lasso10[selDf$id %in% lasso$varList[which(lasso$nMetabDf$lambda==lamb)][[1]]]<-1

# 15 Metabolite:
selDf$lasso15<-0
lamb<-min(lasso$nMetabDf$lambda[!is.na(lasso$nMetabDf$metabs) & lasso$nMetabDf$metabs==15])
selDf$lasso15[selDf$id %in% lasso$varList[which(lasso$nMetabDf$lambda==lamb)][[1]]]<-1

############ Model building ############
# As variable selection:
resDf<-expand.grid(nMetab=c(3,5,7,10,15),varSelect=c("WoACImp","rFImp","lasso"),
                   stringsAsFactors=FALSE)
resDf$ceRF<-resDf$ceGLM<-resDf$ceMultinom<-resDf$misRF<-resDf$misGLM<-resDf$misMultinom<-NA

for(j in 1:nrow(resDf)){
  if(resDf$varSelect[j]!="lasso"){
    vars<-selDf[match(1:resDf$nMetab[j],selDf[,match(resDf$varSelect[j],names(selDf))]),]$id
  }else{
    vars<-selDf$id[as.logical(selDf[,grepl("lasso",names(selDf)) & 
                        resDf$nMetab[j]==gsub("lasso","",names(selDf))])]
  }
  form1<-as.formula(paste("group",paste(vars,collapse="+"),sep="~"))
  
  cv1<-cvFolds(nrow(metab),K=nrow(metab))
  cv2<-data.frame(samp=cv1$subsets,fold=cv1$which)
  
  mis2<-c()
  for(i in 1:nrow(cv2)){
    train<-cv2$samp[cv2$fold!=i]
    test<-cv2$samp[cv2$fold==i]
    
    # Multinomial logit:
    modMult<-multinom(form1,data=metab[train,],
                      MaxNWts=10000,maxit=40,reltol=.0001,trace=FALSE)
    predMult<-log(predict(modMult,newdata=metab[test,],type="probs"))
    predMult[is.infinite(predMult)]<-log(1e-20)
    mm1<-model.matrix(~-1+metab[test,]$group)
    mis2$misMultinom<-c(mis2$misMultinom,
                        as.numeric(!predict(modMult,newdata=metab[test,],"class")==metab$group[test]))
    mis2$ceMultinom<-c(mis2$ceMultinom,-1/nrow(mm1)*sum(mm1*predMult))
    
    # Elastic net:
    tGLMcv<-cv.glmnet(x=as.matrix(metab[train,vars]),alpha=.9,
                      y=metab[train,"group"],family="multinomial",lambda=10**(-seq(0.2,6,.005)),
                      type.measure="deviance",nfolds=10)
    tGLM<-glmnet(x=as.matrix(metab[train,vars]),y=metab[train,"group"],alpha=.9,
                 family="multinomial",lambda=tGLMcv$lambda.min)
    predGLM<-log(predict(tGLM,newx=as.matrix(metab[test,vars]),type="response")[,,1])
    predGLM[is.infinite(predGLM)]<-log(1e-20)
    mis2$misGLM<-c(mis2$misGLM,
                   as.numeric(!levels(metab$group)[which.max(predict(tGLM,newx=as.matrix(metab[test,vars],type="class")))]==
                                metab$group[test]))
    mis2$ceGLM<-c(mis2$ceGLM,-1/nrow(mm1)*sum(mm1*predGLM))
    
    # Random forest:
    modRF<-randomForest(form1,data=metab[train,],ntree=1000)
    predRF<-log(predict(modRF,newdata=metab[test,],type="prob"))
    predRF[is.infinite(predRF)]<-log(1e-20)
    mis2$misRF<-c(mis2$misRF,
                  as.numeric(!predict(modRF,newdata=metab[test,],"class")==metab$group[test]))
    mis2$ceRF<-c(mis2$ceRF,-1/nrow(mm1)*sum(mm1*predRF))
  }
  mis3<-lapply(mis2,mean)
  resDf[j,names(mis3)]<-unlist(mis3)
  print(j)
}
write.csv(resDf,file="resDf.csv",row.names=FALSE)
