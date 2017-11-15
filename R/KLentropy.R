
#' KLentropy
#'
#' Calculates the (weighted) Kozachenko--Leonenko entropy estimator studied in Berrett, Samworth and Yuan (2017), which is based on the \eqn{k}-nearest neighbour distances of the sample.
#'
#' @param x The \eqn{n \times d} data matrix.
#' @param k The tuning parameter that gives the maximum number of neighbours that will be considered by the estimator.
#' @param weights Specifies whether a weighted or unweighted estimator is used. If a weighted estimator is to be used then the default (\code{weights=TRUE}) results in the weights being calculated by \code{\link{L2OptW}}, otherwise the user may specify their own weights.
#'
#' @return The first element of the list is the unweighted estimator for the value of 1 up to the user-specified \eqn{k}. The second element of the list is the weighted estimator, obtained by taking the inner product between the first element of the list and the weight vector.
#'
#' @references \insertRef{BSY2017}{IndepTest}
#' 
#' @examples
#' n=1000; x=rnorm(n); KLentropy(x,30)  # The true value is 0.5*log(2*pi*exp(1)) = 1.42.
#' n=5000; x=matrix(rnorm(4*n),ncol=4)  # The true value is 2*log(2*pi*exp(1)) = 5.68
#' KLentropy(x,30,weights=FALSE)            # Unweighted estimator
#' KLentropy(x,30,weights=TRUE)                           # Weights chosen by L2OptW
#' w=runif(30); w=w/sum(w); KLentropy(x,30,weights=w)  # User-specified weights 
#' 
#' @export
KLentropy=function(x,k,weights=FALSE){
  dim=dim(as.matrix(x)); n=dim[1]; d=dim[2]
  V=pi^(d/2)/gamma(1+d/2)
  if(length(weights)>1){
    w=weights
  }else if(weights==TRUE){
      w=L2OptW(k,d)
  }else{
    w=c(rep(0,k-1),1)
  }
  rho=knn.dist(x,k=k)
  H=(1/n)*colSums(t(log(t(rho)^d*V*(n-1))-digamma(1:k)))
  value=list(); value[[1]]=H; value[[2]]=sum(H*w)
  return(value)
}








