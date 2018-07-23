#' @export
plot.smoothSpl <- function(obj, by = 1 , n = 10, index = NULL, ...){
  xx <- seq(obj$Xcp[1],tail(obj$Xcp,n=1),length.out = obj$NumPoints)
  n <- min(n,dim(obj$Y)[1])
  cols <- rainbow(min(n,30))
  if(is.null(index)) {
    whitch <- seq(1,n,by=by)
  } else {
    whitch <- index
  }
  # Plotting in the clr space fitted curves
  plot.default(xx,obj$Y_clr[1,],ylim=c(min(obj$Y_clr[whitch,]),max(obj$Y_clr[whitch,])),
               type = "n", xlab = "", ylab = "")
  title("Smoothing splines in clr-transformed space")
  for(i in whitch){
    lines(xx,obj$Y_clr[i,], col = cols[i%%length(cols)])
  }

  # Plotting densities in orginal space
  
  plot.default(xx,obj$Y[1,],ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])),
               type = "n", xlab = "", ylab = "")
  title("Density")
  for(i in whitch){
    lines(xx,obj$Y[i,], col = cols[i%%length(cols)])
  }
}
