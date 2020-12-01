
####maize#####
library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

#####DPR######
geno.file<-'pop'
pheno.file<-'pop.fam'
out.file<-'gdd'
Y<-read.table('/mnt/ilustre/centos7users/xiaolong.he/others/GWAS_mengluo/beagle/gdd/data/pop.fam',header =F)
mn<-20;n=2279
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,456,20)
t1<-NULL
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    mat<-cbind(Y[,c(1:5)],y)
    write.table(mat,file = 'pop.fam',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./DPR -b %s -dpr 2 -maf 0.00000001 -w 10000 -s 10000 -o %s",
                   geno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./DPR -b %s -epm output/gdd.param.txt -emu output/gdd.log.txt -predict -o %s",
                   geno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/gdd.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing]
  }
  t1[ij-1] = f(time1)
  writeMat(paste('gdd',ij-1,'_gdd_DPR.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)

}
writeMat(paste('gdd_DPR.mat'),acc1=acc1,acc2=acc2,ytest=ytestc,t1=t1)

#####BSLMM######
geno.file<-'pop'
pheno.file<-'pop.fam'
out.file<-'gdd'
Y<-read.table('/mnt/ilustre/centos7users/xiaolong.he/others/GWAS_mengluo/beagle/gdd/data/pop.fam',header =F)
mn<-20;n=2279
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
t1<-NULL
ytest<-matrix(0,456,20);test<-matrix(0,456,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    mat<-cbind(Y[,c(1:5)],y)
    write.table(mat,file = 'pop.fam',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -b %s -gk 1 -maf 0.00000001 -o %s",
                   geno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -b %s -bslmm 2 -maf 0.00000001 -w 10000 -s 10000 -o %s",
                   geno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -b %s -epm output/gdd.param.txt -emu output/gdd.log.txt -ebv output/gdd.bv.txt -k output/gdd.cXX.txt -predict 1 -o %s",
                   geno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/gdd.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  t1[ij-1] = f(time1)
  writeMat(paste('gdd',ij-1,'_maize_BS.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,test=test,t1=t1)
}
writeMat(paste('gd_BS.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,test=test,t1=t1)



#####BayesR######

geno.file<-'pop'
pheno.file<-'pop.fam'
out.file<-'gdd'
Y<-read.table('/mnt/ilustre/centos7users/xiaolong.he/others/GWAS_mengluo/beagle/gdd/data/pop.fam',header =F)
mn<-20;n=2279
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,456,20)
t1<-NULL
time1 <- Sys.time()
for (ij in 12:21) {
  luoy<-Y[,c(1,6)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    mat<-cbind(Y[,c(1:5)],y)
    write.table(mat,file = 'pop.fam',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./bayesR -bfile %s -out %s ",
                   geno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./bayesR -bfile %s -out %s -predict -model  gdd.model -freq gdd.frq -param gdd.param",
                   geno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='gdd.gv',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing]
  }
  t1[ij-1] = f(time1)
  writeMat(paste('gdd',ij-1,'_maize_BR.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)

}
writeMat(paste('gdd_BR.mat'),acc1=acc1,acc2=acc2,ytest=ytestc,t1=t1)




#####BayesA######
library(R.matlab)
library(glmnet)
library(rrBLUP)
library(BGLR)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

load("mat.RData")
Y<-read.table('pop.fam',header =F)
mn<-20
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
n <- dim(x)[1]
ytest<-matrix(0,length(sample(n,round(n/5),replace=F)),mn)
nIter=10000; burnIn=1000
t1<-NULL
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBA=BGLR(y=ytr,ETA=list( list(X=T,model='BayesA')),
              nIter=nIter,burnIn=burnIn,saveAt='ba')
    acc1[ij-1]<-cor(t%*%fmBA$ETA[[1]]$b+fmBA$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBA$ETA[[1]]$b+fmBA$mu,ytr)
    ytest[,ij-1]<-t%*%fmBA$ETA[[1]]$b+fmBA$mu
    t1[ij-1] = f(time1)
    writeMat(paste('gdd',ij-1,'_BA.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)
  }
}
writeMat(paste('gdd_BA.mat'),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)



#####BayesB######
load("mat.RData")
Y<-read.table('pop.fam',header =F)
mn<-20
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
n <- dim(x)[1]
ytest<-matrix(0,length(sample(n,round(n/5),replace=F)),mn)
nIter=10000; burnIn=1000
t1<-NULL
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBB=BGLR(y=ytr,ETA=list( list(X=T,model='BayesB')),
              nIter=nIter,burnIn=burnIn,saveAt='bb')
    acc1[ij-1]<-cor(t%*%fmBB$ETA[[1]]$b+fmBB$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBB$ETA[[1]]$b+fmBB$mu,ytr)
    ytest[,ij-1]<-t%*%fmBB$ETA[[1]]$b+fmBB$mu
    t1[ij-1] = f(time1)
    writeMat(paste('gdd',ij-1,'_BB.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)
  }
}
writeMat(paste('gdd_BB.mat'),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)



#####Bayes LASSO######
load("mat.RData")
Y<-read.table('pop.fam',header =F)
mn<-20
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
n <- dim(x)[1]
ytest<-matrix(0,length(sample(n,round(n/5),replace=F)),mn)
nIter=10000; burnIn=1000
t1<-NULL
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    fmBL=BGLR(y=ytr,ETA=list( list(X=T,model='BL')),
              nIter=nIter,burnIn=burnIn,saveAt='bl')
    acc1[ij-1]<-cor(t%*%fmBL$ETA[[1]]$b+fmBL$mu,yte)
    acc2[ij-1]<-cor(T%*%fmBL$ETA[[1]]$b+fmBL$mu,ytr)
    ytest[,ij-1]<-t%*%fmBL$ETA[[1]]$b+fmBL$mu
    t1[ij-1] = f(time1)
    writeMat(paste('gdd',ij-1,'_BL.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)
  }
}
writeMat(paste('gdd_BL.mat'),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)




#####rrBLUP######

load("mat.RData")
Y<-read.table('pop.fam',header =F)
mn<-20
acc1<-matrix(0,mn,1)
acc2<-matrix(0,mn,1)
n <- dim(x)[1]
ytest<-matrix(0,length(sample(n,round(n/5),replace=F)),mn)
nIter=10000; burnIn=1000
t1<-NULL
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,6)]
  for (i in 1:1) {
    #luoy<-luoy[sample(nrow(luoy)),]
    #folds <- sample(cut(seq(1,nrow(luoy)),breaks=10,labels=FALSE))
    #ans1<-NULL
    #X<-luoGD[,-1]
    #X<-as.matrix(mat$x)
    M<-as.matrix(x)
    n <- dim(M)[1]
    testing <- sample(n,round(n/5),replace=F)
    trainning <- -testing
    yte <- luoy[testing,2]
    ytr <- luoy[trainning,2]
    #Y<-yte[,2]
    T<-as.matrix(M[trainning,])
    t<-as.matrix(M[testing,])
    RR<- mixed.solve(y=ytr,Z=T)
    acc1[ij-1]<-cor(t%*%RR$u,yte)
    acc2[ij-1]<-cor(T%*%RR$u,ytr)
    ytest[,ij-1]<-t%*%RR$u
    t1[ij-1] = f(time1)
    writeMat(paste('gdd',ij-1,'_RR.mat',sep=""),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)
  }
}
writeMat(paste('gdd_RR.mat'),acc1=acc1,acc2=acc2,ytest=ytest,t1=t1)



