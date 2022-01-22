#'@title oddRatioFunc function
#'
#'
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
#'
#'
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
  a1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(1,1)
  b1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(1,0)
  c1<-CondProb(D,y,z1)$condP+slack
  y[c(i,j)]<-c(0,1)
  d1<-CondProb(D,y,z1)$condP+slack

  return( a1*b1-(c1*d1) )
}

