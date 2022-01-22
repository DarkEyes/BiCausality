#'@title bSCMDepndentGraphFunc function
#'
#'
#'
#'@export
#'
bSCMDepndentGraphFastFunc<-function(mat,IndpThs=0.05)
{
  n<-dim(mat)[1]
  d<-dim(mat)[2]
  # == Create the bootstrapping sequence of D
  D<-VecAlignment(mat)
  #Check dependency of all pairwises
  E0<-matrix(0,nrow=d,ncol=d) # save mean
  E0raw<-E0
  depInfo<-list()

  for(i in seq(d-1))
    for(j in seq(i+1,d))
    {
      str<-sprintf("%d,%d",i,j)
      #print(str)
      bIndpDist<-indpFunc(D,i,j,z=c())

      #check whether i is dependent with j
      E0raw[i,j]<-bIndpDist
      E0raw[j,i]<-bIndpDist
      if(bIndpDist>=IndpThs)
      {
        E0[i,j]<-1
        E0[j,i]<-1
      }
    }
  return(list(E0=E0,E0raw=E0raw) )
}
