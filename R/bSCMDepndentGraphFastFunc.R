#'@title bSCMDepndentGraphFastFunc function
#' @description
#' This function infers dependencies for all pairs of variables without bootstrapping.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#'
#' @return This function returns results of dependency inference among variables.
#' \item{E0}{An adjacency matrix of undirected graph where there is an edge between any pair of variables if they are dependent.}
#' \item{E0raw}{A matrix of the degree of dependency of variable pairs.}

#' @examples
#' bSCMDepndentGraphFastFunc(mat)
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
