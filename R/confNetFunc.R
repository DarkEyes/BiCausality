#'@title confNetFunc function
#' @description
#' This function Computes a confidence network in data mining.
#' Given a set of n transactions or samples in \code{mat} s.t. each transaction has d binary items.
#' The  \code{conf(mat[,j]=1|mat[,i]=1)} is a ratio of a number of samples in jth and ith dimensions that have values equal to
#' one divided by a number of samples in the ith dimension that has a value equal to one.
#' The confNetFunc computes the network where the nodes are dimensions and the edge weights are \code{conf(mat[,j]=1|mat[,i]=1)} for any directed edge from i to j.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param ths is a threshold parameter for cutting of the edge weights. There exists the directed edge from i to j if its edge weight if above or equal \code{ths}.
#' @return This function returns  a binary adjacency matrix \code{confNet} and the weighted adjacency matrix \code{confValMat}.
#' \item{confNet}{A binary adjacency matrix that has \code{confNet[i,j]=1} if \code{confValMat[i,j]>=ths}. Otherwise, it is zero.}
#' \item{confValMat}{A weighted adjacency matrix where \code{confValMat[i,j]} is \code{conf(mat[,j]=1|mat[,i]=1)}.}
#'
#' @examples
#' res<-confNetFunc(mat)
#'
#'@export
#'
confNetFunc<-function(mat,ths=0.1)
{
  d<-dim(mat)[2]
  confNet<-matrix(0,nrow=d,ncol=d)
  D<-VecAlignment(mat)
  confValMat<-matrix(0,nrow=d,ncol=d)
  #===
  for(i in seq(d-1))
    for(j in seq(i+1,d))
    {
      z<-numeric(d)-1
      y<-numeric(d)-1
      y[j]<-c(1)
      z[i]<-c(1)
      a1<-CondProb(D,y=y,z=z)$condP # conf(y|z)
      b1<-CondProb(D,y=z,z=y)$condP # conf(z|y)
      confValMat[i,j]<-a1
      confValMat[j,i]<-b1
      if((a1-b1)>0 && a1>=ths)
        confNet[i,j]<-1
      else if((b1-a1)>0 && b1>= ths)
        confNet[j,i]<-1
      else
        confNet[i,j]<-0
    }
  return(list(confNet=confNet,confValMat=confValMat) )
}
