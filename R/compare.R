xx <- seq(axcp[1],axcp[12],length.out = 100)

myplot <- function(sol, data){
  cumsol <- cumsum(sol[[2]][1,])
  cumdata <- cumsum(data[1,]/100)
  plot(xx,cumsol, type = "l")
  points(axcp,cumdata)
}

myplot(sol,adata)

####################
cut2num <- function(f){ 
  labs <- levels(f) 
  d <- data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), 
                  upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) )) 
  d$midpoints <- rowMeans(d) 
  d 
} 

newdata <- rnorm(1000, 2, 1)
newhist <- cut(newdata, 12, labels = FALSE)
int <- cut2num(cut(newdata, 12))$midpoints
s <- as.vector(table(newhist))

knots <- seq(range(newdata)[1], range(newdata)[2], length.out = 6)
