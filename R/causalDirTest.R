#'@title causalDirTest function
#'
#'
#'
#'@export
#'
causalDirTest<-function(mat,i,j,vi=1,vj=1,z=c(),alpha=0.05,ths = 0.25,nboot=100,slack=0.001)
{
  if(is.null(z))
    z<-numeric(d)-1

  n<-dim(mat)[1]
  bDx<-matrix(0,nboot,n)
  bCausalDirDist<-numeric(nboot)
  for(k in seq(nboot))
  {
    bDx[k,]<-sample(1:n,length(1:n),replace = TRUE)
    nMat<-mat[bDx[k,],]
    D<-VecAlignment(nMat)

    z1<-numeric(d)-1
    y1<-numeric(d)-1
    y1[j]<-c(vj)
    z1[i]<-c(vi)
    a1<-CondProb(D,y=y1,z=z1)$condP
    b1<-CondProb(D,y=z1,z=y1)$condP
    bCausalDirDist[k]<-a1-b1
  }
  testRes<-wilcox.test(x=abs(bCausalDirDist), mu = ths, alternative = "greater")
  confInv<-quantile(bCausalDirDist, c(0+alpha/2, 1-alpha/2))
  bmean <- mean(bCausalDirDist)
  return(list(testRes=testRes,confInv=confInv,bmean=bmean))
}
