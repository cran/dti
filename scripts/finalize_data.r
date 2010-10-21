library(dti)

cat("\nloading data I\n")
bvec <- read.table("/Home/stat/tabelow/DATA/dti_cornell/dvd4/hv-e002093-hardiasset/prepared_data/b-directions.txt")

z <- create.designmatrix.dti(bvec, 1000)

load("/Home/stat/tabelow/DATA/dti_cornell/dvd4/hv-e002093-hardiasset/prepared_data/S0.rsc")
s0 <- final

ttt <- array(0,dim=c(dim(s0),dim(z)[1]))

cat("\nloading data II\n")
for (i in 1:dim(z)[1]) {
  cat(".")
  load(paste("/Home/stat/tabelow/DATA/dti_cornell/dvd4/hv-e002093-hardiasset/prepared_data/S",i,".rsc",sep=""))
  cat(".")
  ttt[,,,i] <- -log(final/s0,base=exp(1))
}

ttt[is.na(ttt)] <- 0
ttt[(ttt==Inf)] <- 0
ttt[(ttt==-Inf)] <- 0

dttt<- dim(ttt)

dim(ttt) <- c(prod(dttt[-4]),dttt[4])
btb<-array(0,c(6,dim(bvec)[1]))
btb[1,]<-bvec[,1]*bvec[,1]
btb[4,]<-bvec[,2]*bvec[,2]
btb[6,]<-bvec[,3]*bvec[,3]
btb[2,]<-2*bvec[,1]*bvec[,2]
btb[3,]<-2*bvec[,1]*bvec[,3]
btb[5,]<-2*bvec[,2]*bvec[,3]
BBB <- btb%*%t(btb)

theta <- solve(BBB)%*% btb %*% t(ttt) / 1000
dim(theta)<-c(6,256,256,72)
eigen <- determine.eigenvalue(aperm(theta,c(2:4,1)))
aniso <- anisotropy(eigen$lambda)

con <- file("unsmoothed_data/Dxx.img","wb")
writeBin(as.vector(theta[1,,,]),con,4)
close(con)
con <- file("unsmoothed_data/Dxy.img","wb")
writeBin(as.vector(theta[2,,,]),con,4)
close(con)
con <- file("unsmoothed_data/Dxz.img","wb")
writeBin(as.vector(theta[3,,,]),con,4)
close(con)
con <- file("unsmoothed_data/Dyy.img","wb")
writeBin(as.vector(theta[4,,,]),con,4)
close(con)
con <- file("unsmoothed_data/Dyz.img","wb")
writeBin(as.vector(theta[5,,,]),con,4)
close(con)
con <- file("unsmoothed_data/Dzz.img","wb")
writeBin(as.vector(theta[6,,,]),con,4)
close(con)

con <- file("unsmoothed_data/D11.img","wb")
writeBin(as.vector(eigen$lambda[,,,3]),con,4)
close(con)
con <- file("unsmoothed_data/D22.img","wb")
writeBin(as.vector(eigen$lambda[,,,2]),con,4)
close(con)
con <- file("unsmoothed_data/D33.img","wb")
writeBin(as.vector(eigen$lambda[,,,1]),con,4)
close(con)

con <- file("unsmoothed_data/Dav.img","wb")
writeBin(as.vector(aniso$trace),con,4)
close(con)
aniso$fa[aniso$fa<0 | aniso$fa>1] <- 0
con <- file("unsmoothed_data/FA.img","wb")
writeBin(as.vector(aniso$fa),con,4)
close(con)

con <- file("unsmoothed_data/EigX.img","wb")
writeBin(as.vector(eigen$theta[,,,3,1]),con,4)
close(con)
con <- file("unsmoothed_data/EigY.img","wb")
writeBin(as.vector(eigen$theta[,,,3,2]),con,4)
close(con)
con <- file("unsmoothed_data/EigZ.img","wb")
writeBin(as.vector(eigen$theta[,,,3,3]),con,4)
close(con)



##### now smoothed results

load("cornell-hardihat5.rsc")

dthat5$theta <- dthat5$theta / 1000
dhat5eigen <- determine.eigenvalue(aperm(dthat5$theta,c(2:4,1)))
dhat5aniso <- anisotropy(dhat5eigen$lambda)

con <- file("smoothed_data/Dxx.img","wb")
writeBin(as.vector(dthat5$theta[1,,,]),con,4)
close(con)
con <- file("smoothed_data/Dxy.img","wb")
writeBin(as.vector(dthat5$theta[2,,,]),con,4)
close(con)
con <- file("smoothed_data/Dxz.img","wb")
writeBin(as.vector(dthat5$theta[3,,,]),con,4)
close(con)
con <- file("smoothed_data/Dyy.img","wb")
writeBin(as.vector(dthat5$theta[4,,,]),con,4)
close(con)
con <- file("smoothed_data/Dyz.img","wb")
writeBin(as.vector(dthat5$theta[5,,,]),con,4)
close(con)
con <- file("smoothed_data/Dzz.img","wb")
writeBin(as.vector(dthat5$theta[6,,,]),con,4)
close(con)

con <- file("smoothed_data/D11.img","wb")
writeBin(as.vector(dhat5eigen$lambda[,,,3]),con,4)
close(con)
con <- file("smoothed_data/D22.img","wb")
writeBin(as.vector(dhat5eigen$lambda[,,,2]),con,4)
close(con)
con <- file("smoothed_data/D33.img","wb")
writeBin(as.vector(dhat5eigen$lambda[,,,1]),con,4)
close(con)

con <- file("smoothed_data/Dav.img","wb")
writeBin(as.vector(dhat5aniso$trace),con,4)
close(con)
dhat5aniso$fa[is.na(dhat5aniso$fa)] <- 0
con <- file("smoothed_data/FA.img","wb")
writeBin(as.vector(dhat5aniso$fa),con,4)
close(con)

con <- file("smoothed_data/EigX.img","wb")
writeBin(as.vector(dthat5$andir[1,,,]),con,4)
close(con)
con <- file("smoothed_data/EigY.img","wb")
writeBin(as.vector(dthat5$andir[2,,,]),con,4)
close(con)
con <- file("smoothed_data/EigZ.img","wb")
writeBin(as.vector(dthat5$andir[3,,,]),con,4)
close(con)
