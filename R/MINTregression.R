

#' MINTregression
#'
#' Performs a goodness-of-fit test of a linear model by testing whether the errors are independent of the covariates.
#'
#' @param x The \eqn{n \times p} design matrix.
#' @param y The response vector of length \eqn{n}.
#' @param k The value of \eqn{k} to be used for estimation of the joint entropy \eqn{H(X,\epsilon)}.
#' @param keps The value of \eqn{k} to be used for estimation of the marginal entropy \eqn{H(\epsilon)}.
#' @param w The weight vector to be used for estimation of the joint entropy \eqn{H(X,\epsilon)}, with the same options as for the \code{\link{KLentropy}} function.
#'@param eps A vector of null errors which should have the same distribution as the errors are assumed to have in the linear model.
#' 
#' @return The \eqn{p}-value corresponding the independence test carried out.
#'
#' @examples
#' \donttest{
#' # Correctly specified linear model
#' x=runif(100,min=-1.5,max=1.5); y=x+rnorm(100)
#' plot(lm(y~x),which=1) 
#' MINTregression(x,y,5,10,w=FALSE,rnorm(10000))
#' # Misspecified mean linear model
#' x=runif(100,min=-1.5,max=1.5); y=x^3+rnorm(100)
#' plot(lm(y~x),which=1)
#' MINTregression(x,y,5,10,w=FALSE,rnorm(10000))
#' # Heteroscedastic linear model
#' x=runif(100,min=-1.5,max=1.5); y=x+x*rnorm(100);
#' plot(lm(y~x),which=1) 
#' MINTregression(x,y,5,10,w=FALSE,rnorm(10000))
#' # Multivariate misspecified mean linear model
#' x=matrix(runif(1500,min=-1.5,max=1.5),ncol=3)
#' y=x[,1]^3+0.3*x[,2]-0.3*x[,3]+rnorm(500)
#' plot(lm(y~x),which=1)
#' MINTregression(x,y,30,50,w=TRUE,rnorm(50000))  
#' }
#' 
#' @references \insertRef{2017arXiv171106642B}{IndepTest} 
#' 
#' @export
MINTregression=function(x,y,k,keps,w=FALSE,eps){
  dim=dim(cbind(x,y)); d=dim[2]; n=dim[1]
  P=x%*%solve(t(x)%*%x)%*%t(x) # Projection matrix  
  B=length(eps)%/%n         # The number of null statistics we can calculate
  nullstat=rep(0,B)
  for(b in 1:B){
    data=eps[((b-1)*n+1):(b*n)]
    res=data-P%*%data; sigma2=sum(res^2)/n
    stdres=res/sqrt(sigma2)
    nullstat[b]=KLentropy(stdres,keps,weights=FALSE)[[2]]-KLentropy(cbind(x,stdres),k,w)[[2]]
  }
  res=y-P%*%y; sigma2=sum(res^2)/n # Residuals and \hat{\sigma}^2
  stdres=res/sqrt(sigma2)
  teststat=KLentropy(stdres,keps,weights=FALSE)[[2]]-KLentropy(cbind(x,stdres),k,w)[[2]]
  p=(1+sum(nullstat>=teststat))/(B+1)
  return(p)
}






