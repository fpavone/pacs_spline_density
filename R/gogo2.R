dyn.load("libpacs.so")
library(Rcpp)
library(RcppEigen)
library(tseries)

cut2num <- function(f){ 
  labs <- levels(f) 
  d <- data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), 
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) )) 
  d$midpoints <- rowMeans(d) 
  d 
} 

ak <- 3
al <- 2
aalpha <- 100
aknots_given <- 1
#adata <- matrix(c(0.3,0.285714286,0.648648648,0.934362935,1.583011583,2.084942085,4.416988413,13.11583012,36.03474903,34.53667954,6.05791506,0),
#                  1,12)
newdata <- rnorm(1000, 2, 1)
newhist <- cut(newdata, 12, labels = FALSE)
s <- t(as.matrix(table(newhist)))/10
int <- cut2num(cut(newdata, 12))$midpoints
#knots <- seq(range(newdata)[1], range(newdata)[2], length.out = 6)
knots <- sort(rnorm(6,2,1))
knots[1] <- range(newdata)[1]
knots[6] <- range(newdata)[2]
adata <- s
axcp <- int
aknots <- knots

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

plot(xx,sol$Y_clr, type = "l")
title(main = aalpha)
points(axcp,sol$Numbers)

save.image(file="normal.RData")
