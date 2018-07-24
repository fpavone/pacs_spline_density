library(splineDensity)
data(particle)
ak <- 3
al <- 2
aalpha <- 10
aknots_given <- 1
axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
u <- log(0.001)
v <- log(200)
classes <- c(u,log(axcp))
midx <- (classes[-1] + classes[-13])/2
lenx <- (classes[-1] - classes[-13])/2
midy <- adata/lenx
aknots <- seq(midx[1],midx[12], length = 7)

test <- function(c,r)
{
  print(c)
  print(r)
  index <- sample(c(1:nrow(midy)),r, replace=T)
  dataset <- midy[index,]/100
  sol <- smoothSplines(ak,al,aalpha,dataset,midx,aknots,cores=c,fast=1)
  sol <- smoothSplines(ak,al,aalpha,dataset,midx,aknots,cores=c,fast=1)
  sol <- smoothSplines(ak,al,aalpha,dataset,midx,aknots,cores=c,fast=1)
  sol <- smoothSplines(ak,al,aalpha,dataset,midx,aknots,cores=c,fast=1)
}

testCV <- function(c,r)
{
  print(c)
  print(r)
  index <- sample(c(1:nrow(midy)),r, replace=T)
  dataset <- midy[index,]/100
  sol <- smoothSplinesVal(ak,al,10^(-6:5),dataset,midx,aknots,cores=c)
  sol <- smoothSplinesVal(ak,al,10^(-6:5),dataset,midx,aknots,cores=c)
  sol <- smoothSplinesVal(ak,al,10^(-6:5),dataset,midx,aknots,cores=c)
  sol <- smoothSplinesVal(ak,al,10^(-6:5),dataset,midx,aknots,cores=c)
}


r <- 1e+06
rCV <- 1e+04
core <- (1:12)

for(i in core){
test(i,r)
testCV(i,rCV)
}
