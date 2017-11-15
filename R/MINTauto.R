
#' MINTauto
#'
#'  Performs an independence test without knowledge of either marginal distribution using permutations and using a data-driven choice of \eqn{k}.
#'
#' @param x The \eqn{n \times d_{X}} data matrix of the \eqn{X} values.
#' @param y The response vector of length \eqn{n \times d_{Y}} data matrix of the \eqn{Y} values.
#' @param kmax The maximum value of \eqn{k} to be considered for estimation of the joint entropy \eqn{H(X,Y)}.
#' @param B1 The number of repetitions used when choosing \eqn{k}, set to 1000 by default.
#' @param B2 The  number of permutations to use for the final test, set at 1000 by default.
#' 
#' @return The \eqn{p}-value corresponding the independence test carried out and the value of \eqn{k} used.
#'
#' @examples
#' \donttest{
#' # Independent univariate normal data
#' x=rnorm(1000); y=rnorm(1000);
#' MINTauto(x,y,kmax=200,B1=100,B2=100)
#' # Dependent univariate normal data
#' library(mvtnorm)
#' data=rmvnorm(1000,sigma=matrix(c(1,0.5,0.5,1),ncol=2))  
#' MINTauto(data[,1],data[,2],kmax=200,B1=100,B2=100)
#' # Dependent multivariate normal data
#' Sigma=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0.5,0,0,0.5,1),ncol=4)
#' data=rmvnorm(1000,sigma=Sigma)
#' MINTauto(data[,1:3],data[,4],kmax=50,B1=100,B2=100)
#'}
#'
#' @references \insertRef{BS2017}{IndepTest}
#' 
#' @export
MINTauto=function(x,y,kmax,B1=1000,B2=1000){
  data=cbind(x,y); d=dim(data)[2]; n=dim(data)[1]
  Hp=matrix(rep(0,B1*kmax),ncol=kmax)
  for(b in 1:B1){
    y1=as.matrix(y)[sample(n),]; data1=cbind(x,y1)
    y2=as.matrix(y)[sample(n),]; data2=cbind(x,y2)
    Hp[b,]=KLentropy(data1,k=kmax)[[1]]-KLentropy(data2,k=kmax)[[1]]
  }
  MSE=(1/B1)*colSums(Hp^2); kopt=which.min(MSE)
  return(c(MINTperm(x,y,kopt,B=B2),kopt))
}

