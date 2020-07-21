#'@title indpFunc function
#'
#'
#'
#'@export
#'
indpFunc<-function(D,i,j,z=c())
{
  if(is.null(z))
    z<-numeric(d)-1

  res<-CondProb(D,y=numeric(d)-1,z=z)
  D<-res$nD
  n<-res$countTotal
  L<-length(D)

  indMagitude<-0

  z1<-numeric(d)-1

  for(i1 in c(0,1) )
    for(j1 in c(0,1))
    {
      y1<-numeric(d) -1
      y1[c(i,j)] <- c(i1,j1)
      y2<-numeric(d) -1
      y2[c(i)] <- i1
      y3<-numeric(d) -1
      y3[c(j)] <- j1

      res2<-CondProb(D,y1,z1)
      n2<-res2$count
      condPair<-res2$condP
      condi<-CondProb(D,y2,z1)$condP
      condj<-CondProb(D,y3,z1)$condP
      indMagitude<-indMagitude+ (abs(condPair - condi*condj)*(n2/n) )
    }

  return(indMagitude)
}
