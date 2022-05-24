#'@title indpFunc function
#'
#' @description
#' This function provides association signs (positive/negative association) inference between i and j.
#' If there is a positive association, it implies i and j trend to have a similar value.
#' For a negative association, however, i and j trend to have an opposite value.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param i is an ith dimension in \code{mat}.
#' @param j is an jth dimension in \code{mat}.
#' @param z is a conditioning d-dimensional vector on \code{mat}.
#' Given k non-negative-bit positions of \code{z}, all k bit positions of samples in the subset of \code{mat} must have similar values with these bits.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction. The default is 0.5.
#' @param nboot is a number of bootstrap replicates for bootstrapping deployed to infer confidence intervals and distributions for hypothesis tests. The default is 100.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#'
#' @return This function returns results of inference of association signs (positive/negative association) between i and j.
#' \item{bmean}{A mean of sign dependency degrees between variables i and j.}
#' \item{confInv}{An \code{alpha}*100th percentile confidence interval of sign dependency degrees between variables i and j.}
#' \item{testRes}{A Mann-Whitney hypothesis test result for an independence test between variables i and j. The null hypothesis is that the distributions of dependency degrees of i,j differ by a location shift of \code{IndpThs} and the alternative is  that distributions of dependency degrees of i,j is shifted greater than \code{IndpThs}. }
#'
#' @examples
#' assocSignTest(mat=mat,i=1,j=2)
#'@export
#'
assocSignTest<-function(mat,i,j,z=c(),alpha=0.05,IndpThs = 0.05,nboot=100)
{
  d<-dim(mat)[2]

  if(is.null(z))
    z<-numeric(d)-1

  n<-dim(mat)[1]
  bDx<-matrix(0,nboot,n)
  bSignDist<-numeric(nboot)
  for(k in seq(nboot))
  {
    bDx[k,]<-sample(1:n,length(1:n),replace = TRUE)
    nMat<-mat[bDx[k,],]
    D<-VecAlignment(nMat)
    bSignDist[k]<-oddDiffFunc(D,i,j,z=z)
  }
  testRes<-wilcox.test(x=abs(bSignDist), mu = IndpThs, alternative = "greater")
  confInv<-quantile(bSignDist, c(0+alpha/2, 1-alpha/2))
  bmean <- mean(bSignDist)
  return(list(testRes=testRes,confInv=confInv,bmean=bmean))
}
