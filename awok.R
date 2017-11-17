############ Prereqs ############
args<-commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=FALSE)
library(methods)
library(nnet)
library(cvTools)
library(doParallel)
library(entropy)
library(dplyr)
library(tidyr)

oldPar<-par()

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

metab0<-metab
pheno0<-pheno
# boot<-sample(1:nrow(metab),replace=TRUE)
# randSub<-sample(1:ncol(metab),size=floor(ncol(metab)/3))
# metab<-metab[boot,randSub]
# pheno<-pheno[boot,]

############ Population class ############
cl<-makeCluster(2)
registerDoParallel(cl)

setClass("population",slots=list(varInclude="list",age="integer",cost="numeric"))
cost<-function(x)
{
  return(unlist(lapply(x,FUN=function(y)
  {
    if(length(table(y))<2)
    {
      return(Inf)
    }
    else
    {
      metab2<-cbind(pheno=pheno$group,as.data.frame(metab[,y]))
      metab2$pheno<-factor(metab2$pheno)
      cv1<-cvFolds(nrow(metab2),K=10,R=10)
      mis<-foreach(j=1:10,.combine=c,.packages=c("nnet"),.inorder=FALSE) %dopar%
      {
        mis2<-c()
        cv2<-data.frame(samp=cv1$subsets[,j],fold=cv1$which)
        for(i in 1:10)
        {
          train<-cv2$samp[cv2$fold!=i]
          test<-cv2$samp[cv2$fold==i]
          
          pred1<-log(predict(multinom(pheno~.,data=metab2[train,],
                                      MaxNWts=10000,maxit=40,reltol=.0001,trace=FALSE),
                             newdata=metab2[test,],type="probs"))
          pred1[is.infinite(pred1)]<-log(1e-200)
          mm1<-model.matrix(~-1+metab2[test,]$pheno)
          mis2<-c(mis2,-1/nrow(mm1)*sum(mm1*pred1))
        }
        mean(mis2)
      }
      return(mean(mis))
    }
  })))
}
#End always run

setGeneric("minCost",function(x,n=1) standardGeneric("minCost"))
setMethod("minCost","population",
          function(x,n)
          {
            return(order(x@cost)[1:n])
          }
)

setGeneric("relFitness",function(x,inverse=FALSE) standardGeneric("relFitness"))
setMethod("relFitness","population",
          function(x,inverse)
          {
            fits<-x@cost
            ecdf1<-ecdf(fits)(fits)
            if(!inverse)
            {
              return(1-ecdf1+min(ecdf1))
            }
            else
            {
              return(ecdf1)
            }
          }
)

setGeneric("recombination",function(x,beta=3) standardGeneric("recombination"))
setMethod("recombination","population",
          function(x,beta)
          {
            relFit<-relFitness(x)
            tinderProb<-qbeta(relFit,shape1=1,shape2=beta)
            tinder<-which(sapply(tinderProb,FUN=function(y) as.logical(rbinom(n=1,size=1,prob=y))))
            weds<-matrix(NA,nrow=0,ncol=2)
            while(length(tinder)>1)
            {
              pair<-sample(tinder,size=2)
              weds<-rbind(weds,pair)
              tinder<-tinder[!tinder %in% pair]
            }
            if(nrow(weds)>1)
            {
              children<-list()
              for(i in 1:nrow(weds))
              {
                mom<-x@varInclude[[weds[i,][1]]]
                dad<-x@varInclude[[weds[i,][2]]]
                # Two random crossover ends:
                cross<-sample(1:(length(mom)-1),2)
                cross<-cross[order(cross)]
                children<-c(children,
                            list(c(dad[1:cross[1]],mom[(cross[1]+1):cross[2]],dad[(cross[2]+1):length(dad)])))
              }
              x@varInclude<-c(x@varInclude,children)
              x@age<-c(x@age,integer(length(children)))
              x@cost<-c(x@cost,cost(children))
            }
            return(x)
          }
)

setGeneric("mutation",function(x) standardGeneric("mutation"))
setMethod("mutation","population",
          function(x)
          {
            relFit<-relFitness(x,inverse=TRUE)
            invFit<-(relFit-min(relFit))/10
            for(j in 1:length(invFit))
            {
              phi<-invFit[j]
              piT<-5/ncol(metab)
              piF<-1-piT
              aTT<-piT+piF*(1-phi)
              aTF<-1-aTT
              aFT<-phi-aTF
              aFF<-1-aFT
              A<-matrix(c(aTT,aTF,aFT,aFF),byrow=TRUE,nrow=2)
              a<-t(A)
              rownames(a)<-colnames(a)<-c("TRUE","FALSE")
              
              trues<-which(x@varInclude[[j]])
              TSwitch<-sample(trues,size=min(rbinom(1,size=length(trues),a['FALSE','TRUE']),
                                             length(trues)))
              x@varInclude[[j]][TSwitch]<-FALSE
              
              falses<-which(!x@varInclude[[j]])
              FSwitch<-sample(falses,size=min(rbinom(1,size=length(falses),a['TRUE','FALSE']),
                                              length(falses)))
              x@varInclude[[j]][FSwitch]<-TRUE
            }
            x@cost<-cost(x@varInclude)
            return(x)
          }
)

setGeneric("death",function(x,popSize) standardGeneric("death"))
setMethod("death","population",
          function(x,popSize)
          {
            costs<-x@cost
            infs<-which(is.infinite(costs))
            if(length(infs)>0)
            {
              x@varInclude<-x@varInclude[-infs]
              x@age<-x@age[-infs]
              x@cost<-x@cost[-infs]
            }
            
            relFit<-relFitness(x)
            ageFit<-ecdf(x@age)(x@age)
            totalFit<-3*relFit+ageFit
            killN<-length(x@varInclude)-popSize
            if(killN>0)
            {
              x@varInclude<-x@varInclude[-order(totalFit)[1:killN]]
              x@age<-x@age[-order(totalFit)[1:killN]]
              x@cost<-x@cost[-order(totalFit)[1:killN]]
            }
            return(x)
          }
)

setGeneric("migration",function(x,varProp,migrationP=0.1) standardGeneric("migration"))
setMethod("migration","population",
          function(x,varProp,migrationP)
          {
            migrants<-ceiling(length(x@varInclude)*migrationP)
            migrants<-birth(migrants,varProp)
            x@varInclude<-c(x@varInclude,migrants@varInclude)
            x@age<-c(x@age,migrants@age)
            x@cost<-c(x@cost,migrants@cost)
            return(x)
          }
)

birth<-function(births,varProp)
{
  success<-rbinom(births,size=ncol(metab),prob=varProp)
  varInclude<-lapply(success,FUN=function(x)
  {
    y<-rep(FALSE,ncol(metab))
    y[sample(1:length(y),x)]<-TRUE
    return(y)
  })
  return(new("population",varInclude=varInclude,age=integer(length(varInclude)),
             cost=cost(varInclude)))
}

generation<-function(popSize=250,varProp=5/ncol(metab),migrationP=.1,iterations=200)
{
  progress<-c()
  births<-popSize
  pop0<-birth(births,varProp)
  progress<-rbind(progress,data.frame(action="start",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
  for(j in 1:iterations)
  {
    pop0@age<-pop0@age+as.integer(1)
    pop0<-recombination(pop0)
    progress<-rbind(progress,data.frame(action="recombination",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    pop0<-mutation(pop0)
    progress<-rbind(progress,data.frame(action="mutation",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    pop0<-migration(pop0,varProp,migrationP=migrationP)
    pop0<-death(pop0,popSize)
    progress<-rbind(progress,data.frame(action="migration",size=length(pop0@varInclude),which=minCost(pop0),val=min(pop0@cost)))
    cat("j is:",j,"\n")
  }
  list(pop=pop0,vars=names(metab),trace=progress)
}

########## Generate Crowds ##########
out<-list()
ptm<-proc.time()
for(i in 1:1)
{
  cat("i is:",i,"\n")
  out[[i]]<-generation()
}
print(proc.time()-ptm)
stopCluster(cl)

assign(paste0("out",args[1]),out)
save(list=paste0("out",args[1]),file=paste0("out",args[1],".RData"))

########### Generate Wisdom ############
library(car)
library(glmnet)
library(ggplot2)
library(gplots)
library(corrplot)
library(colorRamps)

setwd("~/gdrive/Classes/CECS545/FinalProject/")
load("outs.RData")

weights<-c()
for(i in 1:length(outs)) #100
{
  nams<-outs[[i]]$vars
  experts<-outs[[i]]$pop@varInclude[order(outs[[i]]$pop@cost)[1]]
  experts<-sapply(experts,as.numeric)
  rownames(experts)<-nams
  blanks<-matrix(0,nrow=ncol(metab),ncol=1)
  rownames(blanks)<-colnames(metab)
  blanks[match(rownames(experts),rownames(blanks)),]<-experts
  weights<-cbind(weights,apply(blanks,1,mean))
}
key<-key %>% arrange(as.numeric(gsub("M","",id)))
weights<-apply(weights,1,function(x) mean(x))
key$weights<-0
key$weights[match(names(weights),key$id)]<-weights

#Histogram of selection probability:
p1<-ggplot(key %>% filter(weights>0),aes(x=weights,..density..))
#png(file="inclusionP.png",height=3,width=4,units="in",res=600)
p1+geom_histogram(colour="black",fill="grey90",binwidth=.01)+theme_bw()+
  xlab("Inclusion Proportion")+ylab("Density")
#dev.off()

bestKey<-key[key$weights>.025,]

# Abundance distributions:
metab2<-as.matrix(metab[,bestKey$id])
rownames(metab2)<-paste(pheno$group2[match(rownames(metab2),pheno$ptid)],
                        pheno$ptid2[match(rownames(metab2),pheno$ptid)])
colnames(metab2)<-bestKey$biochemical
#png(file="hclust.png",height=6,width=7,units="in",res=600)
plot(hclust(dist(scale(metab2)),method="ward.D"),ann=FALSE)
mtext("Height",2,line=2.25)
#dev.off()

source('~/gdrive/Classes/CECS545/FinalProject/corrplot2.R')

#png(file="corrplot666.png",height=8,width=8,units="in",res=600)
cp<-corrplot2(cor(metab2),order="AOE",tl.cex=.65)
#dev.off()

############ Heatmap / Matrix plot ############
hm1<-heatmap(metab2[,match(rownames(cp),colnames(metab2))],scale="none")
newRows<-c(rownames(metab2)[hm1$rowInd][grepl("sCAD",rownames(metab2)[hm1$rowInd])],
rownames(metab2)[hm1$rowInd][grepl("Non",rownames(metab2)[hm1$rowInd])],
rownames(metab2)[hm1$rowInd][!grepl("Non",rownames(metab2)[hm1$rowInd]) & 
                               !grepl("sCAD",rownames(metab2)[hm1$rowInd])])
metabs3<-metab2[match(newRows,rownames(metab2)),hm1$colInd]
heatmap(metabs3,Rowv=NA,scale="none")

png(file="hm1.png",height=4,width = 3,units="in",res=600)
heatmap(metabs3[1:15,],Rowv=NA,Colv=NA,scale="none",col=rev(matlab.like2(256)),labRow=NA,cexCol=.4)
dev.off()
png(file="hm2.png",height=4,width = 3,units="in",res=600)
heatmap(metabs3[16:27,],Rowv=NA,Colv=NA,scale="none",col=rev(matlab.like2(256)),labRow=NA,cexCol=.4)
dev.off()
png(file="hm3.png",height=4,width = 3,units="in",res=600)
heatmap(metabs3[28:38,],Rowv=NA,Colv=NA,scale="none",col=rev(matlab.like2(256)),labRow=NA,cexCol=.4)
dev.off()


# Boxplots of the selected
metab3<-metab2
metab3$Group<-metab2$pheno
levels(metab3$Group)<-c("Stable CAD","Thrombotic MI","Non-Thrombotic MI")
metab3$pheno<-NULL
metab3<-metab3 %>% gather(id,value,M24:M1031)
metab3<-metab3 %>% left_join(key %>% select(id,biochemical),by="id")
metab3$Group<-factor(metab3$Group,levels=levels(metab3$Group)[c(2,3,1)])
ord<-metab3 %>% group_by(biochemical,Group) %>% dplyr::summarize(median=median(value)) %>%
  arrange(Group,median)
metab3$biochemical<-factor(metab3$biochemical,levels=ord$biochemical)

#png(file="stripCharts.png",height=8,width=8,res=600,units="in")
ggplot(metab3,aes(x=biochemical,y=value,fill=Group,color=Group))+
  theme_bw()+coord_flip()+scale_color_manual(values=c("red","blue","black"))+
  geom_boxplot(fill="white",outlier.colour=NA,position=position_dodge(width=0.7))+
  geom_point(aes(color=Group),position=position_jitterdodge(dodge.width=0.7))+
  xlab("Biochemical")+ylab("Relative Abundance")
#dev.off()

############ Analysis of performance ##############
outT<-c()
for(i in 1:length(outs))
{
  outT<-rbind(outT,outs[[i]]$trace$val)
}
meanPath<-apply(outT,2,mean)

set.seed(33)
samp1<-sample(1:length(outs),10)
#png(file="GAPaths.png",height=4,width=5,units="in",res=600)
plot(1:601,outs[[samp1[1]]]$trace$val,type="l",ylim=c(0,1),xaxt="n",
     main="GA Evolution",xlab="Iterations",ylab="Cross-Entropy Loss",
     col="grey80")
axis(1,at=seq(0,601,75),labels=seq(0,200,25))
for(i in 2:10)
{
  points(1:601,outs[[samp1[i]]]$trace$val,type="l",col="grey80")
}
points(1:601,meanPath,type="l",lwd=2,col="darkblue")
#dev.off()

metab2<-cbind(pheno=pheno$group,as.data.frame(metab[,names(metab)%in%bestKey$id]))
metab2$pheno<-factor(metab2$pheno)

ptm<-proc.time()
compDf<-foreach(i=1:1000,.combine=rbind,.packages="glmnet") %dopar%
{
  # Performance Dataframe
  compDf<-expand.grid(alpha=seq(.25,1,.25),data=c("full","reduced"),mis=NA,n=NA)
  for(j in 1:nrow(compDf))
  {
    set.seed(i)
    if(compDf$data[j]=="full")
    {
      tGLM<-cv.glmnet(x=as.matrix(metab0),alpha=compDf$alpha[j],
                      y=metab2[,"pheno"],family="multinomial",type.measure="class")
      compDf$mis[j]<-min(tGLM$cvm)
      compDf$n[j]<-tGLM$nzero[which.min(tGLM$cvm)]
    }
    else
    {
      tGLM<-cv.glmnet(x=as.matrix(metab2[,names(metab2)!="pheno"]),alpha=compDf$alpha[j],
                      y=metab2[,"pheno"],family="multinomial",type.measure="class")
      compDf$mis[j]<-min(tGLM$cvm)
      compDf$n[j]<-tGLM$nzero[which.min(tGLM$cvm)]
    }
  }
  compDf
}
proc.time()-ptm

sum1<-compDf %>% group_by(alpha,data) %>% 
  summarise(Mis=median(mis),misIQR=IQR(mis),N=mean(n),nSD=sd(n))
#write.csv(sum1,file="sum1.csv")

sum1$mis<-(sum1$Mis/.425)*17
sum1$col<-"darkred"
sum1$col[sum1$data=="reduced"]<-"darkblue"
sum2<-aggregate(n~alpha+data,data=compDf,mean)


library(Hmisc)
source('~/gdrive/Classes/CECS545/FinalProject/back.R')
png(file="b2bMis.png",height=4,width=12,res=300,units="in")
par(mfrow=c(1,4),mar=c(5.1,4.1,4.1,0))
back(compDf$mis[compDf$data=="full" & compDf$alpha==1],compDf$mis[compDf$data=="reduced" & compDf$alpha==1],
             probability=TRUE,brks=seq(0,.425,.025),xlab="",axes=FALSE)
points(c(-.15,.15),sum1$mis[sum1$alpha==1],col=sum1$col[sum1$alpha==1],pch=19)
axis(2,at=0:(length(seq(0,.425,.025))-1),labels=c("0","","","",".10","","","",".20",
                            "","","",".30","","","",".40",""),las=1)
mtext("Lasso")
par(mar=c(5.1,0,4.1,0))
back(compDf$mis[compDf$data=="full" & compDf$alpha==.75],compDf$mis[compDf$data=="reduced" & compDf$alpha==.75],
             probability=TRUE,brks=seq(0,.425,.025),xlab="",axes=FALSE)
points(c(-0.15,.15),sum1$mis[sum1$alpha==.75],col=sum1$col[sum1$alpha==.75],pch=19)
mtext(expression(paste(alpha==.75)))
back(compDf$mis[compDf$data=="full" & compDf$alpha==.5],compDf$mis[compDf$data=="reduced" & compDf$alpha==.5],
             probability=TRUE,brks=seq(0,.425,.025),xlab="",axes=FALSE)
points(c(-0.15,.15),sum1$mis[sum1$alpha==.5],col=sum1$col[sum1$alpha==.5],pch=19)
mtext(expression(paste(alpha==.5)))
par(mar=c(5.1,0,4.1,2.1))
back(compDf$mis[compDf$data=="full" & compDf$alpha==.25],compDf$mis[compDf$data=="reduced"& compDf$alpha==.25],
             probability=TRUE,brks=seq(0,.425,.025),xlab="",axes=FALSE)
points(c(-0.15,.15),sum1$mis[sum1$alpha==.25],col=sum1$col[sum1$alpha==.25],pch=19)
axis(4,at=1:length(seq(0,.425,.025)),labels=c("0","","","",".10",rep("",13)),las=1)
mtext(expression(paste(alpha==.25)))
dev.off()

sum1$Mis<-round(sum1$Mis,3)

set.seed(3)
cvGLM<-cv.glmnet(x=as.matrix(metab2[,names(metab2)!="pheno"]),alpha=.25,
                y=metab2[,"pheno"],family="multinomial",type.measure="class")
cvGLM$nzero[which.min(cvGLM$cvm)]
coefGLM<-coef(cvGLM)
coefGLM<-as.data.frame(cbind(as(coefGLM["sCAD"][[1]],"matrix"),
                             as(coefGLM["Type 1 MI"][[1]],"matrix"),
                             as(coefGLM["Type 2 MI"][[1]],"matrix")))
names(coefGLM)<-c("sCAD","Thrombotic","NonThrombotic")
coefGLM<-coefGLM[rowSums(coefGLM)!=0,]
coefGLM[coefGLM==0]<-NA
coefGLM<-exp(coefGLM)
coefGLM$id<-rownames(coefGLM)
coefGLM<-left_join(coefGLM,key %>%select(id,biochemical,HMDB,platform.x,RI,MASS),by="id")
write.csv(coefGLM,file="coefGLM.csv")

set.seed(2)
sensCV<-cvFolds(nrow(pheno),K=10)
sensCVres<-c()
for(i in 1:10)
{
  train<-sensCV$subsets[sensCV$which!=i]
  test<-sensCV$subsets[sensCV$which==i]
  cvGLM<-cv.glmnet(x=as.matrix(metab2[train,names(metab2)!="pheno"]),alpha=.25,
                   y=metab2[train,"pheno"],family="multinomial",type.measure="class")
  sensCVres<-rbind(sensCVres,cbind(pheno[test,],predict=apply(predict.cv.glmnet(cvGLM,
                                                                                newx=as.matrix(metab2[test,names(metab2)!="pheno"])),1,which.max)))
}
sensCVres$predict<-factor(sensCVres$predict)
levels(sensCVres$predict)<-c("sCAD","Type 1 MI","Type 2 MI")
xtabs(~predict+group,data=sensCVres)


