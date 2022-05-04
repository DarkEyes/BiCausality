#'@title CondProb function
#'
#' @description
#' This function computes a confidence value of \code{y} given \code{c}
#' or \code{conf(y|z)} from an aligned list \code{D}.
#' For any y[i],z[j], their values are -1 by default.
#' The function computes the numbers of transactions
#' that satisfy the following conditions.
#'
#' 1) All transactions must have values at any k position equal to z[k]
#' for any z[k] that is not -1.
#' Let \code{count} be the number of these transactions in \code{D}.
#' 2) All transactions must have values at any k position equal to either z[k] or y[k]
#' that is not -1. Let \code{countTotal} be the number of these transactions in \code{D}.
#'
#'
#'
#' @param D is an aligned list of transactions that was converted from any matrix n by d \code{mat} using
#' \code{D<-VecAlignment(mat)} where n is a number of transactions or samples
#' and d is a number of dimensions for each sample.
#' @param  y is a d-dimensional vector.
#' @param  z is a d-dimensional vector.
#'
#' @return This function returns the ratio \code{condP=count/countTotal}, which is the confidence of \code{y} given \code{z}.
#' \item{condP}{Tthe confidence of \code{y} given \code{z} in \code{D}. }
#' \item{nD}{ The subset of \code{D} such that all transactions
#' have values at any position similar to \code{z[k]} when \code{z[k]} is not -1. }
#' \item{count}{ A number of transactions that have values at any position similar
#' to either \code{z[k]} or \code{y[k]} that is not -1. }
#' \item{countTotal}{ A number of transactions in \code{nD} }
#'
#'@examples
#'d=10 # dimensions of example vectors
#'z<-numeric(d)-1
#'y<-numeric(d)-1
#'y[1]<-c(1)
#'z[c(2,3)]<-c(1,1)
#'CondProb(simData$D,y=y,z=z)$condP # conf(inx1 is 1 |inx 2,3 are 1 ) y|z
#'
#' @export
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
