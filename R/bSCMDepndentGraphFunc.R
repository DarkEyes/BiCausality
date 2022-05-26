#'@title bSCMDepndentGraphFunc function
#'
#' @description
#' This function infers dependencies for all pairs of variables with bootstrapping.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction. The default is 0.5.
#' @param nboot is a number of bootstrap replicates for bootstrapping deployed to infer confidence intervals and distributions for hypothesis tests. The default is 100.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#' @param pflag is a flag for printing progress message (TRUE). The default is FALSE (no printing).
#'
#'
#' @return This function returns results of dependency inference among variables.
#' \item{E0}{An adjacency matrix of undirected graph where there is an edge between any pair of variables if they are dependent.}
#' \item{E0pval}{A matrix of p-values from independence test of pairs of variables.}
#' \item{E0mean}{A matrix of means of dependency degrees between variables.}
#' \item{E0lowbound}{A matrix of lower bounds of dependency-degree confidence intervals between variables.}
#' \item{depInfo\code{[['i,j']]}$bmean}{A mean of dependency degrees between variables i and j.}
#' \item{depInfo\code{[['i,j']]}$confInv}{An \code{alpha}*100th percentile confidence interval of dependency degrees between variables i and j.}
#' \item{depInfo\code{[['i,j']]}$testRes}{A Mann-Whitney hypothesis test result for an independence test between variables i and j. The null hypothesis is that the distributions of dependency degrees of i,j differ by a location shift of \code{IndpThs} and the alternative is  that distributions of dependency degrees of i,j is shifted greater than \code{IndpThs}. }
#' \item{depInfo\code{[['i,j']]}$indices}{A pair of indices of i and j in a numeric vector.}
#' \item{Dboot}{A list of \code{D}s (aligned list of transactions) that are generated from sampling with replacement on input samples (\code{mat}) \code{nboot} times. }
#'
#' @examples
#' \donttest{bSCMDepndentGraphFunc(mat, nboot=50)}
#'
#'@export
#'
bSCMDepndentGraphFunc<-function(mat,nboot=100,alpha=0.05,IndpThs=0.05,pflag=FALSE)
{
  n<-dim(mat)[1]
  d<-dim(mat)[2]
  nMatboot<-list()
  Dboot<-list()
  bDx<-matrix(0,nboot,n)
  # == Create the bootstrapping sequence of D
  for(k in seq(nboot))
  {
    if(pflag==TRUE)
      print(sprintf("bIndpTest-boot#%d",k))
    bDx[k,]<-sample(1:n,length(1:n),replace = TRUE)
    nMatboot[[k]]<-mat[bDx[k,],]
    Dboot[[k]]<-VecAlignment(nMatboot[[k]])
  }
  #Check dependency of all pairwises
  E0<-matrix(0,nrow=d,ncol=d) # save mean
  E0pval<-matrix(0,nrow=d,ncol=d)
  E0mean<-matrix(0,nrow=d,ncol=d)
  E0lowbound<-matrix(0,nrow=d,ncol=d)
  depInfo<-list()

  for(i in seq(d-1))
    for(j in seq(i+1,d))
    {
      str<-sprintf("%d,%d",i,j)
      #print(str)
      bIndpDist<-numeric(nboot)
      for(k in seq(nboot))
      {
        D<-Dboot[[k]]
        bIndpDist[k]<-indpFunc(D,i,j,z=c())
      }
      testRes<-wilcox.test(x=bIndpDist, mu = IndpThs, alternative = "greater")
      confInv<-quantile(bIndpDist, c(0+alpha/2, 1-alpha/2))
      bmean <- mean(bIndpDist)

      #save value in matrices
      E0pval[i,j]<-testRes$p.value
      E0pval[j,i]<-E0pval[i,j]

      E0mean[i,j]<-bmean
      E0mean[j,i]<-E0mean[i,j]

      E0lowbound[i,j]<-confInv[1]
      E0lowbound[j,i]<-E0lowbound[i,j]

      #check whether i is dependent with j
      if(testRes$p.value<=alpha)
      {
        E0[i,j]<-1
        E0[j,i]<-E0[i,j]
        depInfo[[str]]$bmean<-bmean
        depInfo[[str]]$confInv<-confInv
        depInfo[[str]]$testRes<-testRes
        depInfo[[str]]$indices<-c(i,j)
      }
    }
  return(list(E0=E0,E0pval=E0pval,E0mean=E0mean,E0lowbound=E0lowbound,depInfo=depInfo,Dboot=Dboot) )
}
