library(R.matlab)
f <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim08qtn500_va.pheno.txt'
out.file<-'BS_cattle_5024sim08qtn500_va'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.8_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim08qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim08qtn500_va.param.txt -emu output/BS_cattle_5024sim08qtn500_va.log.txt -ebv output/BS_cattle_5024sim08qtn500_va.bv.txt -k output/BS_cattle_5024sim08qtn500_va.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim08qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim08/sim08qtn500_va',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim08/sim08qtn500_va_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim05qtn500_va.pheno.txt'
out.file<-'BS_cattle_5024sim05qtn500_va'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.5_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim05qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim05qtn500_va.param.txt -emu output/BS_cattle_5024sim05qtn500_va.log.txt -ebv output/BS_cattle_5024sim05qtn500_va.bv.txt -k output/BS_cattle_5024sim05qtn500_va.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim05qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim05/sim05qtn500_va',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim05/sim05qtn500_va_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1


geno.file<-'cattle_5024.geno.txt.gz'
pheno.file<-'BS_cattle_5024sim02qtn500_va.pheno.txt'
out.file<-'BS_cattle_5024sim02qtn500_va'
#kin.file<-'BS_cattle_5024sim08qtn500.cXX.txt'
Y<-read.csv('cattle_sim0.2_qtn500_va.csv',header = TRUE)
mn<-20;n=5024
acc1<-matrix(0,20,1)
acc2<-matrix(0,20,1)
ytest<-matrix(0,1005,20);test<-matrix(0,1005,20)
time1 <- Sys.time()
for (ij in 2:21) {
  luoy<-Y[,c(1,ij)]
  for (j in 1:1) {
    y<-matrix(NA,nrow(luoy),1)
    testing <- sample(n,round(n/5),replace=F)
    y[-testing]<-luoy[,2][-testing]
    write.table(y,file = 'BS_cattle_5024sim02qtn500_va.pheno.txt',quote = FALSE,col.names = FALSE,row.names = FALSE)
    system(sprintf("./gemma -g %s -p %s -gk 1 -maf 0.00001 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -bslmm 2 -maf 0.00001 -w 3000 -s 3000 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    system(sprintf("./gemma -g %s -p %s -epm output/BS_cattle_5024sim02qtn500_va.param.txt -emu output/BS_cattle_5024sim02qtn500_va.log.txt -ebv output/BS_cattle_5024sim02qtn500_va.bv.txt -k output/BS_cattle_5024sim02qtn500_va.cXX.txt -predict 1 -o %s",
                   geno.file,pheno.file,out.file),ignore.stdout = FALSE)
    pred<-read.table(file='output/BS_cattle_5024sim02qtn500_va.prdt.txt',header = FALSE)
    acc1[ij-1]<-cor(pred$V1[testing],luoy[,2][testing])
    acc2[ij-1]<-cor(pred$V1[-testing],luoy[,2][-testing])
    ytest[,ij-1]<-pred$V1[testing];test[,ij-1]=testing
  }
  writeMat(paste('BSLMM/sim02/sim02qtn500_va',ij-1,'_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
}
writeMat(paste('BSLMM/sim02/sim02qtn500_va_cattle_BS.mat'),acc1=acc1,acc2=acc2,ytest=ytest,test=test)
t1 = f(time1)
t1
