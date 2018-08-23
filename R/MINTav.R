
#' MINTav
#'
#'  Performs an independence test without knowledge of either marginal distribution using permutations and averaging over a range of values of \eqn{k}.
#'
#' @param x The \eqn{n \times d_{X}} data matrix of the \eqn{X} values.
#' @param y The \eqn{n \times d_{Y}} data matrix of the \eqn{Y} values.
#' @param K The vector of values of \eqn{k} to be considered for estimation of the joint entropy \eqn{H(X,Y)}.
#' @param B The number of permutations to use for the test, set at 1000 by default.
#' 
#' @return The \eqn{p}-value corresponding the independence test carried out.
#'
#' @examples
#' \donttest{
#' # Independent univariate normal data
#' x=rnorm(1000); y=rnorm(1000);
#' MINTav(x,y,K=1:200,B=100)
#' # Dependent univariate normal data
#' library(mvtnorm);
#' data=rmvnorm(1000,sigma=matrix(c(1,0.5,0.5,1),ncol=2))  
#' MINTav(data[,1],data[,2],K=1:200,B=100)
#' # Dependent multivariate normal data
#' Sigma=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0.5,0,0,0.5,1),ncol=4);
#' data=rmvnorm(1000,sigma=Sigma)
#' MINTav(data[,1:3],data[,4],K=1:50,B=100)
#'}
#' 
#' @references \insertRef{2017arXiv171106642B}{IndepTest}
#' 
#' @export
MINTav=function(x,y,K,B=1000){
  data=cbind(x,y); n=dim(data)[1]
  H=KLentropy(data,k=max(K),weights=F)[[1]][K]
  Hp=matrix(rep(0,B*length(K)),nrow=B)
  for(b in 1:B){
    yp=as.matrix(y)[sample(n),]; datap=cbind(x,yp)
    Hp[b,]=KLentropy(datap,k=max(K),weights=F)[[1]][K]
  }
  teststat=sum(H); nullstats=apply(Hp,1,sum)
  p=(1+sum(nullstats<=teststat))/(B+1)
  return(p)
}
