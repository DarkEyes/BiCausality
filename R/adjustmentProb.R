#'@title adjustmentProb function
#'@description
#' This function evaluates the P(Y=yflag|do(X=xflag)) given only marginal distributions using parent adjustment method.
#'
#'@param EValHat is an adjacency matrix of weighted directed causal graph where edge weights are P(Y=yflag|X=xflag) or probabilities of effect being 1 given cause being either 1 for positive association or 0 for negative association.
#'@param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#'@param yflag is value set for Y in P(Y=yflag|X=xflag,z) for the adjustment method.
#'@param xflag is value set for X in P(Y=yflag|X=xflag,z) for the adjustment method.
#'
#'@return This function returns an adjacency matrix of weighted directed causal graph where the edge weights are P(Y=yflag|do(X=xflag) ).
#'
#' @examples
#' adjustmentProb(resC$CausalGRes$EValHat,mat)
#'
#'@export
#'
adjustmentProb<-function(EValHat,mat,yflag=1,xflag=1)
{
  #Work with positive association P(Y=yflag|X=xflag,Z)
  d<-dim(EValHat)[1]
  D<-VecAlignment(mat)
  adjEValHat<-matrix(0,nrow = d, ncol=d)

  for(i in seq(d))
    for(j in seq(d))
    {
      if(EValHat[i,j]>0) # P(Xj|Xi)>0 or xi -> xj
      {
        # if xi has no parent
        np<-sum(EValHat[,i])
        if(np==0)
        {
          adjEValHat[i,j]<-EValHat[i,j]
        }
        else # if xi has any parent(s), it needs adjustment
        {
          filter1<-EValHat[,i]>0 # detect parents of xi (Zi -> xi)

          nBits<-2^np
          #print(filter1)
          for(b in seq(nBits)) # loop to all combination of Zi
          {
            # Z is a parent of xi
            b1<- b-1
            vec1<-num2Bits(b1,n=np)
            z1<-numeric(d)-1
            z1[filter1]<-vec1
            pZi<-CondProb(D,z1,z=numeric(d)-1)$condP #  P(Z=b1)
            message(sprintf("P(Z=%d)=%.2f",b1,pZi))

            y1<-numeric(d)-1
            y1[j]<-c(yflag)
            z2<-z1
            z2[i]<-xflag # set Xi=xflag
            #print(z2)
            # P(Xj=yflag|Xi=xflag,Z=b1)*P(Z=b1)
            a1<-CondProb(D,y=y1,z=z2)$condP
            a1<-ifelse(is.na(a1),0,a1)
            message(sprintf("P(X%d =%d|X%d = %d,z%d)=%.2f",j,yflag,i,xflag,b1,a1))
            adjEValHat[i,j]<-adjEValHat[i,j]+a1*pZi
          }
        }
      }
    }
  return(adjEValHat)
}
#'@title num2Bits function
#'@description
#'Given a natural number and number of bits, the function provides an n-dimensional vector of bits that represents \code{num}.
#' The ith bits of binary vector represents the ith bit of \code{num}.
#' For example, if \code{vec<-num2Bits(num=2,n=4)}, the first bit \code{vec[1]} is 0 and the second bit \code{vec[2]} is 1.
#'
#'@param num is a natural number.
#'@param n is a number of bits representing \code{num}.
#'
#'@return This function returns an n-dimensional vector of bits that represents \code{num}.
#'
#' @examples
#' num2Bits(num=10,n=4)
#'
#'@export
#'
num2Bits<-function(num,n=32)
{
  res <- sapply(num,function(x){ as.integer(intToBits(x))})
  return(as.numeric(res[1:n]) )
}
