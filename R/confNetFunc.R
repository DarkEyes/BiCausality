#'@title confNetFunc function
#' Computing a confidence network in data mining
#' where  there is a weight of edge from i to j as conf(j|i)
#'
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
