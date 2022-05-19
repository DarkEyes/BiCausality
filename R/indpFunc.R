#'@title indpFunc function
#'
#' @description This function computes the degree of dependency between variables.
#' Let i and j be variables, if they are independent, then |p(i,j) -p(i)*p(j)| should be zero.
#' Given the samples in the n by d matrix \code{mat} where n is a number of samples and d is a number of dimensions,
#' an aligned list of transactions \code{D} is computed by
#' \code{D<-VecAlignment(mat)}.
#'
#'@param D is an aligned list of transactions that was converted from \code{mat}.
#'@param i is an ith dimension in \code{mat}.
#'@param j is an jth dimension in \code{mat}.
#'
#'@return This function returns the degree of dependency between variables:
#' zero implies both variables are independent, and non-zero value implies the degree of dependency (higher implies more dependent degree).
#'
#' @examples
#' indpFunc(simData$D,i=1,j=2)
#'
#'@export
#'
indpFunc<-function(D,i,j,z=c())
{
  d<-length(D[[1]]$name)
  if(is.null(z)) # go full D without conditional variables
    z<-numeric(d)-1

  res<-CondProb(D,y=numeric(d)-1,z=z) #get the total n
  D<-res$nD
  n<-res$countTotal
  L<-length(D)

  indMagitude<-0

  z1<-numeric(d)-1

  for(i1 in c(0,1) )
    for(j1 in c(0,1))
    {
      y1<-numeric(d) -1
      y1[c(i,j)] <- c(i1,j1) # supp(i,j)
      y2<-numeric(d) -1
      y2[c(i)] <- i1 #supp(i)
      y3<-numeric(d) -1
      y3[c(j)] <- j1 #supp(j)

      res2<-CondProb(D,y1,z1)
      n2<-res2$count
      condPair<-res2$condP
      condi<-CondProb(D,y2,z1)$condP
      condj<-CondProb(D,y3,z1)$condP
      # |p(i,j) -p(i)*p(j)|*weight
      indMagitude<-indMagitude+ (abs(condPair - condi*condj)*(n2/n) )
    }

  return(indMagitude)
}
