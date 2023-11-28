#'@title CausalGraphInferMainFunc function
#' @description
#' A framework to infer causality on binary data using techniques in frequent pattern mining and estimation statistics. Given a set of individual vectors S=\{x\} where x(i) is a realization value of binary variable i, the framework infers empirical causal relations of binary variables i,j from S in a form of causal graph G=(V,E) where V is a set of nodes representing binary variables and there is an edge from i to j in E if the variable i causes j. The framework determines dependency among variables as well as analyzing confounding factors before deciding whether i causes j.
#'
#' Note that all statistics (e.g. means) and confidence intervals as well as hypothesis testing are inferred by bootstrapping.
#'
#' @param mat is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction. The default is 0.5.
#' @param nboot is a number of bootstrap replicates for bootstrapping deployed to infer confidence intervals and distributions for hypothesis tests. The default is 100.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#' @param CausalThs is a threshold for the degree of causal direction In the causal-direction test, to claim that any variables have causal relations, the degree of causal direction must greater than this value significantly. The default is 0.1.
#'
#' @return This function returns causal inference results. #TODO: provide list of results.
#' \item{depRes}{The result of inferring dependencies between all pairs of variables.}
#' \item{ConfoundRes}{The result of filtering associations without true causal directions from any confounding factor.}
#' \item{CausalGRes}{The result of inferring causal directions between all pairs of dependent variables that have no confounding factors.}
#' \item{depRes$E0}{An adjacency matrix of undirected graph where there is an edge between any pair of variables if they are dependent.}
#' \item{depRes$E0pval}{A matrix of p-values from independence test of pairs of variables.}
#' \item{depRes$E0mean}{A matrix of means of dependency degrees between variables.}
#' \item{depRes$E0lowbound}{A matrix of lower bounds of dependency-degree confidence intervals between variables.}
#' \item{depRes$depInfo$'i,j'$bmean}{A mean of dependency degrees between variables i and j.}
#' \item{depRes$depInfo$'i,j'$confInv}{An \code{alpha}*100th percentile confidence interval of dependency degrees between variables i and j.}
#' \item{depRes$depInfo$'i,j'$testRes}{A Mann-Whitney hypothesis test result for an independence test between variables i and j. The null hypothesis is that the distributions of dependency degrees of i,j differ by a location shift of \code{IndpThs} and the alternative is  that distributions of dependency degrees of i,j is shifted greater than \code{IndpThs}. }
#' \item{depRes$depInfo$'i,j'$indices}{A pair of indices of i and j in a numeric vector.}
#' \item{depRes$Dboot}{A list of \code{D}s (aligned list of transactions) that are generated from sampling with replacement on input samples (\code{mat}) \code{nboot} times. }
#' \item{ConfoundRes$E1}{An adjacency matrix of undirected graph after filtering associations without true causal directions from any confounding factor.}
#' \item{ConfoundRes$E2}{A matrix of associations that have confounding factors where \code{E2[i,j]=0} if no confounding factor and \code{E2[i,j]=k} if k is a confounding factor of i and j.}
#' \item{CausalGRes$Ehat}{An adjacency matrix of directed causal graph where \code{CausalGRes$Ehat[i,j]=1} implies i causes j.}
#' \item{CausalGRes$EValHat}{An adjacency matrix of weighted directed causal graph where edge weights are estimated means of probabilities of effect being 1 given cause being either 1 for positive association or 0 for negative association using CondProb() and bootstrapping to estimate}
#' \item{CausalGRes$causalInfo$'i,j'$CDirConfValInv}{An \code{alpha}*100th percentile confidence interval of estimated conditional probability of effect j being 1/0 given cause i's value being either the same (positive association) or opposite (negative association).}
#' \item{CausalGRes$causalInfo$'i,j'$CDirConfInv}{An \code{alpha}*100th percentile confidence interval of estimated causal direction degree of i cause j. }
#' \item{CausalGRes$causalInfo$'i,j'$CDirmean}{A mean-estimated-causal-direction degree of i cause j.}
#' \item{CausalGRes$causalInfo$'i,j'$testRes2}{A Mann-Whitney hypothesis test result for existence of causal direction. The null hypothesis is that the distributions of causal-direction degrees of i,j differ by a location shift of \code{CausalThs} and the alternative is  that distributions of causal-direction degrees of i,j is shifted greater than \code{CausalThs}.}
#' \item{CausalGRes$causalInfo$'i,j'$testRes1}{A Mann-Whitney hypothesis test result for existence of association by odd differences from \code{oddDiffFunc()}. The null hypothesis is that the distributions of absolute odd difference of i,j differ by a location shift of \code{IndpThs} and the alternative is  that distributions of absolute odd difference of i,j is shifted greater than \code{IndpThs}.}
#' \item{CausalGRes$causalInfo$'i,j'$sign}{A direction of i,j association: 1 for positive, 0 for negative, and -1 for no association.}
#' \item{CausalGRes$causalInfo$'i,j'$SignConfInv}{An \code{alpha}*100th percentile confidence interval of i,j odd difference from bootstrapping. }
#' \item{CausalGRes$causalInfo$'i,j'$Signmean}{A mean of i,j odd difference from bootstrapping.}
#'
#' @examples
#' \donttest{resC<-CausalGraphInferMainFunc(mat = mat, nboot =50)}
#'
#'@export
#'
CausalGraphInferMainFunc<-function(mat,alpha=0.05,nboot=100,IndpThs=0.05,CausalThs = 0.1)
{
  message("Inferring dependent graph")
  res<-bSCMDepndentGraphFunc(mat,nboot=nboot,alpha=alpha,IndpThs=IndpThs)
  message("Removing confounder(s)")
  res2<-bSCMdeConfoundingGraphFunc(res,IndpThs=IndpThs,alpha=alpha)
  message("Inferring causal graph")
  res3<-bSCMCausalGraphFunc(res2$E1,res$Dboot ,alpha=alpha,SignThs=IndpThs,CausalThs = CausalThs)
  return(list(depRes = res, ConfoundRes= res2, CausalGRes= res3))
}
