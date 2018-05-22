dyn.load("libpacs.so")
library(Rcpp)
library(RcppEigen)
library(tseries)

## SCRIPT FOR DATA GIVEN BY PROF. MENAFOGLIO

ak <- 3
al <- 2
aalpha <- 10
aknots_given <- 1
adata <- read.matrix("../input/data", sep = " ", skip = 1)
axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
#aknots <- c(0.063,0.25,1,4,16,63,100)


u <- log(0.001)
v <- log(200)
classes <- c(u,log(axcp))
midx <- (classes[-1] + classes[-13])/2
lenx <- (classes[-1] - classes[-13])/2
midy <- adata/lenx
aknots <- seq(midx[1],midx[12], length = 7)

fun <- function(k,l,alpha,data,xcp,knots,num_points = 100) 
{
  # Checking if data is a matrix
  if ( !is.matrix(data) ) 
  {
    err <- simpleError("data must be a matrix type.")
    stop(err)
  }

  # Creating equispaced knots if not given
  if( length(knots) == 1 )
  {
    u <- xcp[1]
    v <- tail(xcp,n=1)
    size <- knots
    step <- (v - u)/(size-1)
    knots_ <- seq(u,v, by = step)
    .Call("mymain",as.integer(k),as.integer(l),alpha,data,xcp,knots_,as.integer(num_points))
  }
  else 
    .Call("mymain",as.integer(k),as.integer(l),alpha,data,xcp,knots,as.integer(num_points))
}

nrow <- dim(adata)[1]


sol <- fun(ak,al,aalpha,as.data.frame(midy/100),midx,aknots)



xx <- seq(midx[1],midx[12],length.out = 100)

myplot <- function(sol, data){
  cumsol <- cumsum(sol$Y[nrow,])/sum(sol$Y[nrow,])
  cumdata <- cumsum(data[nrow,]/100)
  plot(xx,cumsol, type = "l")
  points(midx,cumdata)
}

myplot(sol,adata)
# Plotting in the clr space fitted curves with input data as points
plot(xx,sol$Y_clr[nrow,],ylim=c(min(sol$Numbers),max(sol$Numbers)), type = "l")
title(main = aalpha)
points(midx,sol$Numbers,pch = 19 ,col = "tomato3")

cols <- rainbow(30)
# Plotting in the clr space fitted curves
plot(xx,sol$Y_clr[nrow,],ylim=c(min(sol$Numbers),max(sol$Numbers)), type = "l")
for(i in seq(1,(nrow-1),by=3)){
  lines(xx,sol$Y_clr[i,], col = cols[i%%length(cols)])
}

# Plotting densities in orginal space
plot(xx,sol$Y[nrow,],ylim=c(min(midy[nrow,]),max(midy[nrow,]))/100, type = "l")
for(i in seq(1,(nrow-1),by=3)){
  lines(xx,sol$Y[i,], col = cols[i%%length(cols)])
}

save.image(file="sol.RData")

