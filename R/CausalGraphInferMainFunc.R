#'@title CausalGraphInferMainFunc function
#' @description
#' A framework to infer causality on binary data using techniques in frequent pattern mining and estimation statistics. Given a set of individual vectors S={x} where x(i) is a realization value of binary variable i, the framework infers empirical causal relations of binary variables i,j from S in a form of causal graph G=(V,E) where V is a set of nodes representing binary variables and there is an edge from i to j in E if the variable i causes j. The framework determines dependency among variables as well as analyzing confounding factors before deciding whether i causes j.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction.
#' @param nboot is a number of bootstrap replicates for bootstrapping deployed to infer confidence intervals and distributions for hypothesis tests.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the degree of dependency must greater than this value significantly.
#' @param CausalThs is a threshold for the degree of causal direction In the causal-direction test, to claim that any variables have causal relations, the degree of causal direction must greater than this value significantly.
#'
#' @return This function returns causal inference results. #TODO: provide list of results
#'
#' @examples
#' #resC<-BiCausality::CausalGraphInferMainFunc(mat = simData$mat,CausalThs=0.1, nboot =50, IndpThs=0.05)
#'
#'@export
#'
CausalGraphInferMainFunc<-function(mat,alpha=0.05,nboot=100,IndpThs=0.05,CausalThs = 0.1,slack=0.001)
{
  print("Inferring dependent graph")
  res<-bSCMDepndentGraphFunc(mat,nboot=nboot,alpha=alpha,IndpThs=IndpThs)
  print("Removing confounder(s)")
  res2<-bSCMdeConfoundingGraphFunc(res,IndpThs=IndpThs,alpha=alpha)
  print("Inferring causal graph")
  res3<-bSCMCausalGraphFunc(res2$E1,res$Dboot ,alpha=alpha,SignThs=IndpThs,CausalThs = CausalThs,slack=slack)
  return(list(depRes = res, ConfoundRes= res2, CausalGRes= res3))
}
