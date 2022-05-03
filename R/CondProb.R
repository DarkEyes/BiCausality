#'@title CondProb function
#' The function serves as a confidence function in data mining conf(y|z)
#' Given a set of transaction D
#'
#'
#'@export
#'
CondProb<-function(D,y,z)
{
  p<-0
  filter<-z != -1
  if(sum(filter) == 0) # go full D without conditional variables
    nD<-D
  else # Keep only binary that are match with the given z.
  {
    nD<-list()
    for(i in seq(length(D)) )  #For each unique binary D[[i]]
    {
      flag<-sum(D[[i]]$name[filter] == z[filter]) # check whether D[[i]] bits in z's positions are the same as z's
      if(flag == sum(filter)) # if so, add the D[[i]] into nD
      {
        nD[[names(D)[i] ]] <- D[[i]]
      }
    }
  }
  filterY<- y!= -1
  count<-0
  countTotal<-0
  for(i in seq(length(nD)) ) #For each unique binary nD[[i]] in nD
  {
    flag<-sum(nD[[i]]$name[filterY] == y[filterY])
    countTotal<-countTotal+nD[[i]]$count
    if(flag == sum(filterY))  #  check whether nD[[i]] bits in y's positions are the same as y's
    {
      count<-count+nD[[i]]$count
    }
  }

  return(list("condP"=count/countTotal, nD=nD,countTotal=countTotal,count=count ))
}
