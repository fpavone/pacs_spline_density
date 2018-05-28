dyn.load("libpacs.so")
library(Rcpp)
library(RcppEigen)
library(tseries)
source("../R/mypackage.R")

## SCRIPT FOR DATA GIVEN BY PROF. MENAFOGLIO

ak <- 3
al <- 2
aalpha <- 10
aknots_given <- 1
adata <- read.matrix("../input/data", sep = " ", skip = 1)
axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)

u <- log(0.001)
v <- log(200)
classes <- c(u,log(axcp))
midx <- (classes[-1] + classes[-13])/2
lenx <- (classes[-1] - classes[-13])/2
midy <- adata/lenx
aknots <- seq(midx[1],midx[12], length = 7)


sol <- smoothingSplines(ak,al,aalpha,midy/100,midx,aknots)
plot(sol)



save.image(file="sol.RData")


## CUMULATIVE PLOT
# xx <- seq(midx[1],midx[12],length.out = 100)
#
# myplot <- function(sol, data){
#   cumsol <- cumsum(sol$Y[nrow,])/sum(sol$Y[nrow,])
#   cumdata <- cumsum(data[nrow,]/100)
#   plot(xx,cumsol, type = "l")
#   points(midx,cumdata)
# }
#
# myplot(sol,adata)

## PLOT IN CLR SPACE WITH GIVEN DATA
# # Plotting in the clr space fitted curves with input data as points
# plot(xx,sol$Y_clr[nrow,],ylim=c(min(sol$Numbers),max(sol$Numbers)), type = "l")
# title(main = aalpha)
# points(midx,sol$Numbers,pch = 19 ,col = "tomato3")
