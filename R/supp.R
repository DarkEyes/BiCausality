#'@title supp function
#'
#' @description
#' This function computes a support value from a matrix \code{X} given a \code{values}.
#'
#' @param X is a matrix n by d where n is a number of transactions or samples
#' and d is a number of dimensions for each sample.
#' @param  values is a d-dimensional vector
#' we use to count how many of it within \code{X}.
#'
#' @return This function returns the support of \code{values} in \code{X} by counting
#'  the ratio of how many samples in \code{X} are similar to \code{values}
#'@examples
#' x <- rbinom(n=100, size=1, prob=0.5)
#' ny<-rbinom(n=100, size=1, prob=0.25)
#' y <- x |  ny
#' supp(X=cbind(x,y),values=c(1,1) )
#'
#' @export
#'
supp<-function(X,values)
{
  count<-0
  n<-0
  flag=0
  if(is.null(dim(X)[1]))
  {
    n<-length(X)
    flag=1
  }
  else
    n<-dim(X)[1]
  for(i in seq(n))
  {
    if(flag==1)
      row<-X[i]
    else
      row<-X[i,]
    if(sum(row==values) == length(values))
      count=count+1
  }
  return(count/n)
}
