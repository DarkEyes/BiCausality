#' @title  comparePredAdjMatrix2TrueAdjMat
#'
#' @description
#'
#' comparePredAdjMatrix2TrueAdjMat is a support function that can compare two adjacency matrices: ground-truth and inferred matrices.
#'
#' @param trueAdjMat a ground-truth matrix.
#' @param adjMat an inferred matrix.
#'
#' @return This function returns a list of precision \code{prec}, recall \code{rec}, and F1 score \code{F1} of inferred vs. groundtruth matrices.
#'
#' @examples
#' ## Generate simulation data
#' #G<-matrix(FALSE,10,10) # groundtruth
#' #G[1,c(4,7,8,10)]<-TRUE
#' #G[2,c(5,7,9,10)]<-TRUE
#' #G[3,c(6,8,9,10)]<-TRUE
#' #comparePredAdjMatrix2TrueAdjMat(trueAdjMat=G,adjMat=G)
#'
#'@export
comparePredAdjMatrix2TrueAdjMat<-function(trueAdjMat,adjMat)
{
  TP<-0
  FP<-0
  FN<-0
  for(i in seq(10))
    for(j in seq(10))
    {
      if(trueAdjMat[i,j] && adjMat[i,j])
        TP<-TP+1
      else if( (!trueAdjMat[i,j]) && adjMat[i,j])
        FP<-FP+1
      else if(trueAdjMat[i,j] && (!adjMat[i,j]) )
        FN<-FN+1
    }
  prec<-TP/(TP+FP)
  rec<-TP/(TP+FN)
  F1<-2*prec*rec/(prec+rec)
  return(list(prec=prec,rec=rec,F1=F1))
}
