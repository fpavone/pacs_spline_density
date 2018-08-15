remove.packages("splineDensity")

devtools::install_github("fpavone/pacs_spline_density", ref = "ROpenMP")

library(splineDensity)

example(smoothSplines)

# Add minor tick marks
library(Hmisc)

##########################
##### plot chapter 3 #####
##########################
# particle-size data

# all

whitch<-1:406
cols <- rainbow(40)
plot.default(xx,obj$Y[1,],ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
                type = "n", xlab = "", ylab = "")
plot.default(xx,obj$Y[1,],xlim=c(-3.5,4.2),ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
                type = "n", xlab = "log(Particle-size diameter)          [log(mm)]", ylab = "")
cols[1]
for(i in whitch){
     lines(xx,obj$Y[i,], col = cols[i%%length(cols) + 1])
}

minor.tick(nx=4, tick.ratio=0.75)

# just a few

whitch <- seq(1,12,by=1)
cols <- rainbow(8)
plot.default(xx,obj$Y[1,],ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
             type = "n", xlab = "", ylab = "")
plot.default(xx,obj$Y[1,],xlim=c(-3.5,4.2),ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
             type = "n", xlab = "log(Particle-size diameter)          [log(mm)]", ylab = "")
cols[1]
for(i in whitch){
  lines(xx,obj$Y[i,], col = cols[i%%length(cols) + 1])
}

i<-i+1;lines(xx,obj$Y[i,], col = cols[i%%length(cols) + 1])

# just a meaningful few

whitch <- c(2,3,4,8, 10, 32)
cols <- rainbow(6)
plot.default(xx,obj$Y[1,],xlim=c(-3.5,4.2),ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
             type = "n", xlab = "log(Particle-size diameter)          [log(mm)]", ylab = "")
cols[1]
j<-1
for(i in whitch){
  lines(xx,obj$Y[i,], col = cols[j+1])
  j<-j+1
}
lines(xx,obj$Y[32,], col = "forestgreen")
lines(xx,obj$Y[2,], col = "#FF3300")

# or

plot.default(xx,obj$Y[1,],xlim=c(-3.5,4.2),ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
             type = "n", xlab = "log(Particle-size diameter)          [log(mm)]", ylab = "")
cols[1]
j<-1
for(i in whitch){
  lines(xx,obj$Y[i,], col = cols[j])
  j<-j+1
}




