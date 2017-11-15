 
#' L2OptW
#'
#' Calculates a weight vector to be used for the weighted Kozachenko--Leonenko estimator. The weight vector has minimum \eqn{L_2} norm subject to the linear and sum-to-one constraints of (2) in Berrett, Samworth and Yuan (2017). 
#'
#' @param k The tuning parameter that gives the number of neighbours that will be considered by the weighted Kozachenko--Leonenko estimator.
#' @param d The dimension of the data.
#'
#' @return The weight vector that is the solution of the optimisation problem.
#'
#' @examples
#' # When d < 4 there are no linear constraints and the returned vector is (0,0,...,0,1).
#' L2OptW(100,3)    
#' w=L2OptW(100,4)
#' plot(w,type="l")
#' w=L2OptW(100,8);
#' # For each multiple of 4 that d increases an extra constraint is added.
#' plot(w,type="l")  
#' w=L2OptW(100,12)
#' plot(w, type="l") # This can be seen in the shape of the plot
#'
#' @references \insertRef{BSY2017}{IndepTest}
#' 
#' @export
L2OptW=function(k,d){
  dprime=floor(d/4)
  if(dprime==0){
    return(c(rep(0,k-1),1))
    }else{
    G=matrix(rep(0,(dprime+1)*k),ncol=k)
    G[1,]=rep(1,k)
    for(l in 1:dprime){
      G[l+1,]=exp(lgamma(1:k+2*l/d)-lgamma(1:k)) # exp(lgamma-lgamma) better than gamma/gamma for large k
    }
    A=G%*%t(G)
    Lambda=solve(A,c(1,rep(0,dprime)))
    return(as.vector(t(G)%*%Lambda))
  }
}
