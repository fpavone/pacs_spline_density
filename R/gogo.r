dyn.load("libpacs.so")
library(Rcpp)
library(RcppEigen)
library(tseries)

## SCRIPT FOR DATA GIVEN BY PROF. MENAFOGLIO

ak <- 3
al <- 2
aalpha <- 0.25
aknots_given <- 1
adata <- read.matrix("../input/data", sep = " ", skip = 1)
axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
aknots <- c(0.063,0.25,1,4,16,63,100)

fun <- function(k,l,alpha,knots_given,data,xcp,knots) {
  .Call("mymain",as.integer(k),as.integer(l),alpha,as.integer(knots_given),data,xcp,knots)
}



sol <- fun(ak,al,aalpha,aknots_given,adata,axcp,aknots)



xx <- seq(axcp[1],axcp[12],length.out = 100)

myplot <- function(sol, data){
  cumsol <- cumsum(sol[[2]][1,])
  cumdata <- cumsum(data[1,]/100)
  plot(xx,cumsol, type = "l")
  points(axcp,cumdata)
}

myplot(sol,adata)

plot(xx,sol$Y_clr[dim(adata)[1],], type = "l")
title(main = aalpha)
points(axcp,sol$Numbers)

save.image(file="sol.RData")

