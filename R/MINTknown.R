#' MINTknown
#'
#' Performs an independence test when it is assumed that the marginal distribution of \eqn{Y} is known and can be simulated from.
#'
#' @param x The \eqn{n \times d_X} data matrix of \eqn{X} values.
#' @param y The \eqn{n \times d_Y} data matrix of \eqn{Y} values.
#' @param k The value of \eqn{k} to be used for estimation of the joint entropy \eqn{H(X,Y)}.
#' @param ky The value of \eqn{k} to be used for estimation of the marginal entropy \eqn{H(Y)}.
#' @param w The weight vector to used for estimation of the joint entropy \eqn{H(X,Y)}, with the same options as for the \code{\link{KLentropy}} function.
#' @param wy The weight vector to used for estimation of the marginal entropy \eqn{H(Y)}, with the same options as for the \code{\link{KLentropy}} function.
#' @param y0 The data matrix of simulated \eqn{Y} values.
#' 
#' @return The \eqn{p}-value corresponding the independence test carried out.
#'
#' @examples
#' library(mvtnorm)
#' x=rnorm(1000); y=rnorm(1000);
#' # Independent univariate normal data
#' MINTknown(x,y,k=20,ky=30,y0=rnorm(100000))  
#' library(mvtnorm)
#' # Dependent univariate normal data
#' data=rmvnorm(1000,sigma=matrix(c(1,0.5,0.5,1),ncol=2))
#' # Dependent multivariate normal data
#' MINTknown(data[,1],data[,2],k=20,ky=30,y0=rnorm(100000))   
#' Sigma=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0.5,0,0,0.5,1),ncol=4)
#' data=rmvnorm(1000,sigma=Sigma)
#' MINTknown(data[,1:3],data[,4],k=20,ky=30,w=TRUE,wy=FALSE,y0=rnorm(100000))
#'
#' @references \insertRef{BS2017}{IndepTest}
#' 
#' @export
MINTknown=function(x,y,k,ky,w=FALSE,wy=FALSE,y0){
  n=dim(cbind(x,y))[1]
  B=dim(as.matrix(y0))[1]%/%n   # The number of null statistics we can calculate
  nullstat=rep(0,B)
  for(b in 1:B){
    yb=as.matrix(y0)[((b-1)*n+1):(b*n),] # Extracting the bth null sample
    nullstat[b]=KLentropy(yb,k=ky,weights=wy)[[2]]-KLentropy(cbind(x,yb),k=k,weights=w)[[2]] # Calculating null stat
  }
  data=cbind(x,y)          # Real data
  stat=KLentropy(y,k=ky,weights=wy)[[2]]-KLentropy(data,k=k,weights=w)[[2]] # Calculating the real test stat
  p=(1+sum(nullstat>=stat))/(B+1) # Comparing the real test stat to the null test stats
  return(p)
}



