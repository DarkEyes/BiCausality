#'@title indpFunc function
#'
#'
#'
#'@export
#'
assocSignTest<-function(mat,i,j,z=c(),alpha=0.05,ths = 0.05,nboot=100,slack=0.001)
{
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
    bSignDist[k]<-oddDiffFunc(D,i,j,z=z,slack=slack)
  }
  testRes<-wilcox.test(x=abs(bSignDist), mu = ths, alternative = "greater")
  confInv<-quantile(bSignDist, c(0+alpha/2, 1-alpha/2))
  bmean <- mean(bSignDist)
  return(list(testRes=testRes,confInv=confInv,bmean=bmean))
}
