

#' MINTknown
#'
#' Performs an independence test without knowledge of either marginal distribution using permutations.
#'
#' @param x The \eqn{n \times d_X} data matrix of \eqn{X} values.
#' @param y The \eqn{n \times d_Y} data matrix of \eqn{Y} values.
#' @param k The value of \eqn{k} to be used for estimation of the joint entropy \eqn{H(X,Y)}.
#' @param w The weight vector to used for estimation of the joint entropy \eqn{H(X,Y)}, with the same options as for the \code{\link{KLentropy}} function.
#' @param B The number of permutations to use, set at 1000 by default.
#' 
#' @return The \eqn{p}-value corresponding the independence test carried out.
#'
#' @examples
#' # Independent univariate normal data
#' x=rnorm(1000); y=rnorm(1000)
#' MINTperm(x,y,k=20,B=100)
#' # Dependent univariate normal data
#' library(mvtnorm)
#' data=rmvnorm(1000,sigma=matrix(c(1,0.5,0.5,1),ncol=2))  
#' MINTperm(data[,1],data[,2],k=20,B=100)
#' # Dependent multivariate normal data
#' Sigma=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0.5,0,0,0.5,1),ncol=4)
#' data=rmvnorm(1000,sigma=Sigma)
#' MINTperm(data[,1:3],data[,4],k=20,w=TRUE,B=100)
#'
#' @references \insertRef{BS2017}{IndepTest}
#' 
#' @export
MINTperm=function(x,y,k,w=FALSE,B=1000){
  data=cbind(x,y); n=dim(data)[1]
  Hp=rep(0,B)
  for(b in 1:B){
    yp=as.matrix(y)[sample(n),]
    Hp[b]=KLentropy(cbind(x,yp),k,weights=w)[[2]]
  }
  H=KLentropy(data,k,weights=w)[[2]]
  p=(1+sum(Hp<=H))/(B+1)
  return(p)
}






