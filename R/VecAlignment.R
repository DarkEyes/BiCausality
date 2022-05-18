
#'@title VecAlignment function
#'
#' @description
#' This function rearranges the samples in the \code{mat} into
#' an aligned list of transactions, which is mainly used by other functions in the package.
#' Suppose mat[i,] is a binary vector  we are interested, we use  \code{A<-bin2dec(mat[i,])}
#' to store the decimal value of \code{mat[i,]} in \code{A}. Then, we call \code{D[[A]]$count}
#' to get number of samples in \code{mat} that are similar to \code{mat[i,]} and
#' the \code{D[[A]]$name} is \code{mat[i,]}.
#'
#'
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#'
#' @return This function returns the list \code{D}, is an aligned list of transactions that was converted from any matrix n by d \code{mat}.
#'
#' @examples
#' VecAlignment(mat=simData$mat)
#'
#'@export
#'
VecAlignment<-function(mat)
{
  # mat must have the dimension n by d where n and d must be greater than 1.
  D<-list()
  d<-dim(mat)[1]
  for(i in seq(d))
  {
    D[[sprintf("%d",bin2dec(mat[i,]))]]$name<-mat[i,]
    if(is.null(D[[sprintf("%d",bin2dec(mat[i,]))]]$count)==TRUE )
      D[[sprintf("%d",bin2dec(mat[i,]))]]$count<-1
    else
      D[[sprintf("%d",bin2dec(mat[i,]))]]$count<-newMat[[sprintf("%d",bin2dec(mat[i,]))]]$count+1
  }
  return(D)
}
#'@title bin2dec function
#'
#'@description
#'This function convertes a binary vector into its decimal value.
#'@param x is a binary vector where \code{X[i]} is the ith bit of vector.
#'
#'@return This function returns a decimal value of \code{X}.
#'
#'@examples
#'bin2dec(x=c(1,1,1,0))
#'
#'@export
#'
bin2dec<-function(x)
{
  newx<-0
  d<-length(x)
  for(i in seq(d))
  {
    newx <- newx+x[i]*2^(d-i)
  }
  return(newx)
}
