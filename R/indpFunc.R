#'@title indpFunc function
#'  |p(i,j) -p(i)*p(j)| should be zero if i and j are independent
#'
#'
#'@export
#'
indpFunc<-function(D,i,j,d,z=c())
{
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
