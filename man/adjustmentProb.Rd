% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustmentProb.R
\name{adjustmentProb}
\alias{adjustmentProb}
\title{adjustmentProb function}
\usage{
adjustmentProb(EValHat, mat, yflag = 1, xflag = 1)
}
\arguments{
\item{EValHat}{is an adjacency matrix of weighted directed causal graph where edge weights are P(Y=yflag|X=xflag) or probabilities of effect being 1 given cause being either 1 for positive association or 0 for negative association.}

\item{mat}{is a matrix n by d where n is a number of transactions or samples and d is a number of dimensions.}

\item{yflag}{is value set for Y in P(Y=yflag|X=xflag,z) for the adjustment method.}

\item{xflag}{is value set for X in P(Y=yflag|X=xflag,z) for the adjustment method.}
}
\value{
This function returns an adjacency matrix of weighted directed causal graph where the edge weights are P(Y=yflag|do(X=xflag) ).
}
\description{
This function evaluates the P(Y=yflag|do(X=xflag)) given only marginal distributions using parent adjustment method.
}
\examples{
adjustmentProb(resC$CausalGRes$EValHat,mat)

}
