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
SepalLengthCm <- iris$Sepal.Length
Species <- iris$Species

iris1 <- SepalLengthCm[iris$Species==levels(iris$Species)[1]]
h1 <- hist(iris1, freq = FALSE, nclass = 12, plot = F)

lines(density(iris1))

xx <- seq(ssol$Xcp[1],tail(ssol$Xcp,n=1),length.out = ssol$NumPoints)
lines(xx,ssol$Y[1,], col = 'red')



iris2 <- SepalLengthCm[iris$Species==levels(iris$Species)[2]]
h2 <- hist(iris2, freq = FALSE, nclass = 12)

iris3 <- SepalLengthCm[iris$Species==levels(iris$Species)[3]]
h3 <- hist(iris3, freq = FALSE, nclass = 12)




library(splineDensity)

k <- 3
l <- 2
alpha <- 100
# uu <- 4.25
# vv <- 5.85

mmidx <- h1$mids
mmidy <- matrix(c(h1$density, h1$density),nrow=2, ncol = 15, byrow=T)
knots <- 7 # or as an alternative,  seq(4, 6, length = 7)
ssol <- smoothingSplines(k,l,alpha, mmidy,mmidx,knots)
plot(ssol)


h1 <- hist(iris1, freq = FALSE, nclass = 12)
lines(density(iris1))
xx <- seq(ssol$Xcp[1],tail(ssol$Xcp,n=1),length.out = ssol$NumPoints)
lines(xx,ssol$Y[1,], col = 'red')

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
