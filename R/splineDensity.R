#' Estimate density from histogram
#'
#' @param k smoothing splines degree
#' @param l order of derivative in the penalization term
#' @param alpha weight for penalization
#' @param data an object of class "matrix" containing data to be smoothed, row by row
#' @param xcp vector of control points
#' @param knots either vector of knots for the splines or a integer for the number of equispaced knots
#' @param num_points number of points of the grid where to evaluate the density estimated
#' @param prior prior used for zero-replacements. This must be one of "perks", "jeffreys", "bayes_laplace", "sq" or "default"
#' @param cores number of cores for parallel execution, if the option was enabled before installing the ppackage
#' @param fast 1 if maximal performance is required (print statements suppressed), 0 otherwise
#' @return An object of class "smoothSpl"
#' @description Given raw (discretized) distributional observations, smoothSplines computes the density
#' function that 'best' fits data, in a trade-off between smooth and least squares approximation. 
#' @details The original discretized densities are not directly smoothed. The centred logratio transformation is
#' first applied, as defined
#' @references J. Machalová, K. Hron & G.S. Monti (2016): 
#' Preprocessing of centred logratio transformed density functions 
#' using smoothing splines. Journal of Applied Statistics, 43:8, 1419-1435.
#' @examples
#' library(splineDensity)
#' data(particle)
#' ak <- 3
#' al <- 2
#' aalpha <- 10
#' aknots_given <- 1
#' axcp <- c(0.063,0.125,0.25,0.5,1,2,4,8,16,31.5,63,100)
#' u <- log(0.001)
#' v <- log(200)
#' classes <- c(u,log(axcp))
#' midx <- (classes[-1] + classes[-13])/2
#' lenx <- (classes[-1] - classes[-13])/2
#' midy <- adata/lenx
#' aknots <- seq(midx[1],midx[12], length = 7)
#' sol <- smoothSplines(ak,al,aalpha,midy/100,midx,aknots)
#' @useDynLib splineDensity
#' @export 
#' 

smoothSplines <- function(k,l,alpha,data,xcp,knots,num_points = 100, prior = "default", cores = 1, fast = 0)
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
    obj <- .Call("smoothingSplines_",as.integer(k),as.integer(l),alpha,
                 data,xcp,knots_,as.integer(num_points),as.integer(prior_num), as.integer(cores), as.integer(fast))
  }
  else
    obj <- .Call("smoothingSplines_",as.integer(k),as.integer(l),alpha,
                 data,xcp,knots,as.integer(num_points),as.integer(prior_num), as.integer(cores), as.integer(fast))
  
  class(obj) <- "smoothSpl"
  return(obj)
}

