#setwd("~/Desktop/VM/pacs_spline_density")

clr_transf <- function(y)
{
  a = psych::geometric.mean(y)
  return(log(y/a))
}

data<- as.matrix(read.table("input/data"))
knots<-c(0.063, 0.25, 1, 4, 16, 63, 100)
k = 3+1
l = 2
alpha = 100
xcp <-data[1,]
ycp1 <-c(-0.213133, -0.894851, -0.630159, -0.142864, -0.175535, -0.139299, 0.201232, 0.566227, 1.09655, 1.70559, 0.478097, -1.85185)
ycp2 <-c(-0.223162, -0.567666, -0.114948, 0.637918, 0.261027, 0.405066, 0.828658, 1.2187, 1.41222, 2.02723, -1.5092, -1.5092 )
ycp3 <-c(-0.213854, -0.312588, -0.151722, 0.435619, 0.351849, 0.0211847, 0.255081, 1.35102, 2.22775, 1.97813, -1.36658, -1.36658)
ycp <- c(0.3,0.285714286,0.648648648,0.934362935,1.583011583,2.084942085,4.416988413,13.11583012,36.03474903,34.53667954,6.05791506,0.012)
w = c(1,1,1,1,1,1,1,1,1,1,1,1)
ch=2
load("ssp_clr.Rdata")
#ssp_clr(knots,xcp,ycp1,w,k,l,alpha,ch)
#ssp_clr(knots,xcp,ycp2,w,k,l,alpha,ch)
#ssp_clr(knots,xcp,ycp3,w,k,l,alpha,ch)

teval <- seq(xcp[1], xcp[12], length.out = 102)
x <- knots
t <- xcp
f <- clr_transf(ycp)
der <- l
alfa <- alpha
nknots <- length(knots)

out <- ssp_clr(log(knots),log(xcp),clr_transf(ycp),w,k,l,alpha,ch, log(teval))

