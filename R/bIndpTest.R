#'@title bIndpTest function
#'
#'
#'
#'@export
#'
bIndpTest<-function(mat,i,j,z=c(),alpha=0.05,ths = 0.05,nboot=100,pflag=FALSE)
{
  d<-dim(mat)[2]
  if(is.null(z))
    z<-numeric(d)-1
  #z[2]<- -1

  n<-dim(mat)[1]
  bDx<-matrix(0,nboot,n)
  bIndpDist<-numeric(nboot)
  for(k in seq(nboot))
  {
    if(pflag==TRUE)
      print(sprintf("bIndpTest-boot#%d",k))
    bDx[k,]<-sample(1:n,length(1:n),replace = TRUE) #sampling index vector
    nMat<-mat[bDx[k,],] #samplimng from mat and save to nMat
    D<-VecAlignment(nMat) # align
    bIndpDist[k]<-indpFunc(D,i,j,z=z)
  }
  testRes<-wilcox.test(x=bIndpDist, mu = ths, alternative = "greater")
  confInv<-quantile(bIndpDist, c(0+alpha/2, 1-alpha/2))
  bmean <- mean(bIndpDist)
  return(list(testRes=testRes,confInv=confInv,bmean=bmean))
}
