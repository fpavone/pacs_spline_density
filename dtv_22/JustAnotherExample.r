#---------------------------
######  EXAMPLE #2  ########
#####  IRIS DATASET  #######
#---------------------------

#---------------------------
#     Density estimation
#---------------------------

SepalLengthCm <- iris$Sepal.Length
Species <- iris$Species

# you can do this

# kernel (not kernAl)

ggplot(iris, aes(x=SepalLengthCm, colour=Species, fill=Species)) +
    geom_density(alpha=.3) +
    xlim(3.5,8.5)+
    xlab("Sepal Length (cm)") +
    ylab("Density")


# OR

# The sm.density.compare( ) function in the sm package allows you
# to superimpose the kernal density plots of two or more groups.
library(sm)

# plot densities
sm.density.compare(SepalLengthCm, Species, xlab="Sepal Length (cm)")
title(main="Sepal Length Distribution by Species")

# add legend
colfill<-c(2:(2+length(levels(Species))))
legend(x = "topright", levels(Species), fill=colfill)



#---------------------------
#    Histograms
#---------------------------
ggplot(data=iris, aes(x=SepalLengthCm))+
  geom_histogram(binwidth=0.2, color="black", aes(fill=Species)) +
  xlab("Sepal Length (cm)") +
  ylab("Frequency") +
  ggtitle("Histogram of Sepal Length") +
  theme(plot.title = element_text(hjust = 0.5))

# or avoiding "stratification":

# separate plots
ggplot(data=iris, aes(x=SepalLengthCm)) +
  geom_histogram() +
  facet_grid(~Species)

# same plot
ggplot(data=iris,aes(x=SepalLengthCm,group=Species,fill=Species)) +
  geom_histogram(position="dodge",binwidth=0.25) +
  theme_bw()

# same plot, overlappping
ggplot(data=iris,aes(x=SepalLengthCm,group=Species,fill=Species)) +
  geom_histogram(position="identity",alpha=0.5,binwidth=0.2) +
  theme_bw()


#---------------------------
#    Now using our code
#---------------------------

# preparing histograms - to feed our algorithm
SepalLengthCm <- iris$Sepal.Length
Species <- iris$Species

iris1 <- SepalLengthCm[iris$Species==levels(iris$Species)[1]]
h1 <- hist(iris1, nclass = 12, plot = F)

iris2 <- SepalLengthCm[iris$Species==levels(iris$Species)[2]]
h2 <- hist(iris2, nclass = 12, plot = F)

iris3 <- SepalLengthCm[iris$Species==levels(iris$Species)[3]]
h3 <- hist(iris3, nclass = 12, plot = F)

# setting parameters
library(splineDensity)

k <- 3
l <- 2
alpha <- 100
# uu <- 4.25
# vv <- 5.85

#### Iris: level 1
midx1 <- h1$mids
midy1 <- matrix(h1$density, nrow=1, ncol = length(h1$density), byrow=T)
knots <- 7 # or as an alternative,  seq(4, 6, length = 7)
sol1 <- smoothSplines(k,l,alpha,midy1,midx1,knots)
plot(sol1)

h1 <- hist(iris1, freq = FALSE, nclass = 12)
lines(density(iris1))
xx1 <- seq(sol1$Xcp[1],tail(sol1$Xcp,n=1),length.out = sol1$NumPoints)
lines(xx1,sol1$Y[1,], col = 'red')

#### Iris: level 2
midx2 <- h2$mids
midy2 <- matrix(h2$density, nrow=1, ncol = length(h2$density), byrow=T)
knots <- 7 # or as an alternative,  seq(4, 6, length = 7)
sol2 <- smoothSplines(k,l,alpha,midy2,midx2,knots)
plot(sol2)

h2 <- hist(iris2, freq = FALSE, nclass = 12)
lines(density(iris2))
xx2 <- seq(sol2$Xcp[1],tail(sol2$Xcp,n=1),length.out = sol2$NumPoints)
lines(xx2,sol2$Y[1,], col = 'red')

#### Iris: level 3
midx3 <- h3$mids
midy3 <- matrix(h3$density, nrow=1, ncol = length(h3$density), byrow=T)
knots <- 7 # or as an alternative,  seq(4, 6, length = 7)
sol3 <- smoothSplines(k,l,alpha,midy3,midx3,knots)
plot(sol3)

h3 <- hist(iris3, freq = FALSE, nclass = 12)
lines(density(iris3))
xx3 <- seq(sol3$Xcp[1],tail(sol3$Xcp,n=1),length.out = sol3$NumPoints)
lines(xx3,sol3$Y[1,], col = 'blue')

solCv <- smoothSplinesVal(k,l,10^(-5:5), midy3, midx3, knots)


# just for completeness, the lines of code of the example
# modified for taking into account a prescribed number of rows
# data(particle)
# ak <- 3
# al <- 2
# aalpha <- 10
# axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
# u <- log(0.001)
# v <- log(200)
# classes <- c(u,log(axcp))
# midx <- (classes[-1] + classes[-13])/2
# lenx <- (classes[-1] - classes[-13])/2
# midy <- adata/lenx
# aknots <- seq(midx[1],midx[12], length = 7)
# sol <- smoothSplines(ak,al,aalpha,matrix(midy[1,], nrow=1, ncol=12)/100,midx,aknots)
# plot(sol)
