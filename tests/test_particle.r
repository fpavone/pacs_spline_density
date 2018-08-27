library(splineDensity)

data(particle)

x_cp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)

u <- log(0.001)
v <- log(200)
classes <- c(u,log(x_cp))

midx <- (classes[-1] + classes[-13])/2
lenx <- (classes[-1] - classes[-13])/2
midy <- adata/lenx
aknots <- seq(midx[1],midx[12], length = 7)

sol <- smoothSplines(k=3,l=2,alpha=10,midy/100,midx,aknots)
plot(sol)

print("All is working! :)")