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
    v <- utils::tail(xcp,n=1)
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
