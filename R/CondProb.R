#'@title CondProb function
#'
#'
#'
#'@export
#'
CondProb<-function(D,y,z)
{
  p<-0
  filter<-z != -1
  if(sum(filter) == 0)
    nD<-D
  else
  {
    nD<-list()
    for(i in seq(length(D)) )
    {
      flag<-sum(D[[i]]$name[filter] == z[filter])
      if(flag == sum(filter))
      {
        nD[[names(D)[i] ]] <- D[[i]]
      }
    }
  }
  filterY<- y!= -1
  count<-0
  countTotal<-0
  for(i in seq(length(nD)) )
  {
    flag<-sum(nD[[i]]$name[filterY] == y[filterY])
    countTotal<-countTotal+nD[[i]]$count
    if(flag == sum(filterY))
    {
      count<-count+nD[[i]]$count
    }
  }

  return(list("condP"=count/countTotal, nD=nD,countTotal=countTotal,count=count ))
}
