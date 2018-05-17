dyn.load("libpacs.so")
library(Rcpp)
library(RcppEigen)
library(tseries)

ak <- 3
al <- 2
aalpha <- 1
aknots_given <- 1
#adata <- matrix(c(0.3,0.285714286,0.648648648,0.934362935,1.583011583,2.084942085,4.416988413,13.11583012,36.03474903,34.53667954,6.05791506,0),
#                  1,12)
adata <- read.matrix("../input/data", sep = " ", skip = 1)
axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
aknots <- c(0.063,0.25,1,4,16,63,100)

fun <- function(k,l,alpha,knots_given,data,xcp,knots) {
  .Call("mymain",as.integer(k),as.integer(l),alpha,as.integer(knots_given),data,xcp,knots)
}



sol <- fun(ak,al,aalpha,aknots_given,adata,axcp,aknots)
