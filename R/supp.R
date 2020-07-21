#'@title supp function
#'
#'
#'
#'@export
#'
supp<-function(X,values)
{
  count<-0
  n<-0
  flag=0
  if(is.null(dim(X)[1]))
  {
    n<-length(X)
    flag=1
  }
  else
    n<-dim(X)[1]
  for(i in seq(n))
  {
    if(flag==1)
      row<-X[i]
    else
      row<-X[i,]
    if(sum(row==values) == length(values))
      count=count+1
  }
  return(count/n)
}
