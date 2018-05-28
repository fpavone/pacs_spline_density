smoothingSplines <- function(k,l,alpha,data,xcp,knots,num_points = 100, prior = "default") 
{
  # Checking if data is a matrix
  if ( !is.matrix(data) ) 
  {
    err <- simpleError("data must be a matrix type.")
    stop(err)
  }
  
  # Converting prior to numeric type
  prior_num <- 0
  if ( prior == "perks" ) prior_num <- 1
  else if ( prior == "jeffreys" ) prior_num <- 2
  else if ( prior == "bayes_laplace" ) prior_num <- 3
  else if ( prior == "sq" ) prior_num <- 4

  
  # Creating equispaced knots if not given
  if( length(knots) == 1 )
  {
    u <- xcp[1]
    v <- tail(xcp,n=1)
    size <- knots
    step <- (v - u)/(size-1)
    knots_ <- seq(u,v, by = step)
    obj <- .Call("mymain",as.integer(k),as.integer(l),alpha,
                 data,xcp,knots_,as.integer(num_points),as.integer(prior_num))
  }
  else 
   obj <- .Call("mymain",as.integer(k),as.integer(l),alpha,
                data,xcp,knots,as.integer(num_points),as.integer(prior_num))
  
  class(obj) <- "smoothSpl"
  return(obj)
}

plot.smoothSpl <- function(obj, by = 1 , n = 10,...){
  xx <- seq(obj$Xcp[1],tail(obj$Xcp,n=1),length.out = obj$NumPoints)
  n <- min(n,dim(obj$Y)[1])
  cols <- rainbow(min(n,30))
  whitch <- seq(1,n,by=by)
  # Plotting in the clr space fitted curves
  plot.default(xx,obj$Y_clr[1,],ylim=c(min(obj$Y_clr[whitch,]),max(obj$Y_clr[whitch,])), 
               type = "l", xlab = "", ylab = "")
  title("Smoothing splines in clr-transformed space")
  for(i in whitch){
    lines(xx,obj$Y_clr[i,], col = cols[i%%length(cols)])
  }
  
  # Plotting densities in orginal space
  plot.default(xx,obj$Y[1,],ylim=c(min(obj$Y[whitch,]),max(obj$Y[whitch,])), 
               type = "l", xlab = "", ylab = "")
  title("Density")
  for(i in whitch){
    lines(xx,obj$Y[i,], col = cols[i%%length(cols)])
  }
  
}
