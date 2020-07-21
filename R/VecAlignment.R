
#'@title VecAlignment function
#'
#'
#'
#'@export
#'
VecAlignment<-function(mat)
{
  # mat must have the dimension n by d where n and d must be greater than 1.
  newMat<-list()
  d<-dim(mat)[1]
  for(i in seq(d))
  {
    newMat[[sprintf("%d",bin2dec(mat[i,]))]]$name<-mat[i,]
    if(is.null(newMat[[sprintf("%d",bin2dec(mat[i,]))]]$count)==TRUE )
      newMat[[sprintf("%d",bin2dec(mat[i,]))]]$count<-1
    else
      newMat[[sprintf("%d",bin2dec(mat[i,]))]]$count<-newMat[[sprintf("%d",bin2dec(mat[i,]))]]$count+1
  }
  return(newMat)
}
#'@title bin2dec function
#'
#'
#'
#'@export
#'
bin2dec<-function(x)
{
  newx<-0
  d<-length(x)
  for(i in seq(d))
  {
    newx <- newx+x[i]*2^(d-i)
  }
  return(newx)
}
