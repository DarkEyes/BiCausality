#'@title oddRatioFunc function
#'
#'@description
#'Given the samples in the n by d matrix \code{mat} where n is a number of samples and d is a number of dimensions.
#'This function computes an odd ratio value of variables of ith and jth dimensions from
#'a given an aligned list of transactions \code{D} (compute by \code{D<-VecAlignment(mat)}).
#'
#'@param D is an aligned list of transactions that was converted from \code{mat}.
#'@param i is an ith dimension in \code{mat} for computing the odd ratio with.
#'@param j is an jth dimension in \code{mat} for computing compute the odd ratio with.
#'@param z is a conditioning d-dimensional vector on \code{D}.
#' Given k non-negative-bit positions of \code{z}, all k bit positions of samples in the subset of \code{D} must have similar values with these bits.
#'@param slack is a parameter to prevent the issue of division by zero.
#'
#'@return This function returns an odd ratio value of variables of ith and jth dimensions from \code{D}.
#'
#' @examples
#' oddRatioFunc(simData$D,i=1,j=2)
#'
#'@export
#'
oddRatioFunc<-function(D,i,j,z=c(),slack=0.001)
{
  d<-length(D[[1]]$name)
  if(is.null(z))
    z<-numeric(d)-1

  res<-CondProb(D,y=numeric(d)-1,z=z)
  D<-res$nD
  n<-res$countTotal
  L<-length(D)

  oddMagitude<-0

  z1<-numeric(d)-1
  y<-numeric(d)-1

  y[c(i,j)]<-c(0,0)
  a1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(1,1)
  b1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(1,0)
  c1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(0,1)
  d1<-CondProb(D,y,z1)$condP+slack

  return(a1*b1/(c1*d1))
}

#'@title oddDiffFunc function
#'@description
#'Given the samples in the n by d matrix \code{mat} where n is a number of samples and d is a number of dimensions.
#'This function computes an odd difference value of variables of ith and jth dimensions from
#'a given an aligned list of transactions \code{D} (compute by \code{D<-VecAlignment(mat)}).
#'
#'@param D is an aligned list of transactions that was converted from \code{mat}.
#'@param i is an ith dimension in \code{mat} for computing the odd difference with.
#'@param j is an jth dimension in \code{mat} for computing compute the odd difference with.
#'@param z is a conditioning d-dimensional vector on \code{D}.
#' Given k non-negative-bit positions of \code{z}, all k bit positions of samples in the subset of \code{D} must have similar values with these bits.
#'
#'@return This function returns an odd difference value of variables of ith and jth dimensions from \code{D}.
#'
#' @examples
#' oddDiffFunc(simData$D,i=1,j=2)
#'
#'@export
#'
oddDiffFunc<-function(D,i,j,z=c(),slack=0.001)
{
  d<-length(D[[1]]$name)
  if(is.null(z))
    z<-numeric(d)-1

  res<-CondProb(D,y=numeric(d)-1,z=z)
  D<-res$nD
  n<-res$countTotal
  L<-length(D)

  oddMagitude<-0

  z1<-numeric(d)-1
  y<-numeric(d)-1

  y[c(i,j)]<-c(0,0)
  a1<-CondProb(D,y,z1)$condP
  y[c(i,j)]<-c(1,1)
  b1<-CondProb(D,y,z1)$condP
  y[c(i,j)]<-c(1,0)
  c1<-CondProb(D,y,z1)$condP
  y[c(i,j)]<-c(0,1)
  d1<-CondProb(D,y,z1)$condP

  return( a1*b1-(c1*d1) )
}

