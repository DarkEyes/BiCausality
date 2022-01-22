#'@title bSCMDepndentGraphFunc function
#'
#'
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
