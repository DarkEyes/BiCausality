#'@title bSCMCausalGraphFunc function
#'
#' @description
#' This function infers a causal graph from a result of confounding factor filtering by \code{bSCMdeConfoundingGraphFunc()}.
#'
#' @param E1 is an adjacency matrix of undirected graph after filtering associations without true causal directions from any confounding factor.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction. The default is 0.5.
#' @param Dboot is a list of \code{D}s (aligned list of transactions) that are generated from sampling with replacement on input samples (\code{mat}) \code{nboot} times.
#' @param SignThs is a threshold for the degree of dependency for association direction inference. In the independence test of sign direction, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#' @param CausalThs is a threshold for the degree of causal direction In the causal-direction test, to claim that any variables have causal relations, the degree of causal direction must greater than this value significantly. The default is 0.1.
#'
#' @return This function returns causal inference results from E1 matrix that is an output of \code{bSCMdeConfoundingGraphFunc}.
#' \item{Ehat}{An adjacency matrix of directed causal graph where \code{CausalGRes$Ehat[i,j]=1} implies i causes j.}
#' \item{EValHat}{An adjacency matrix of weighted directed causal graph where edge weights are estimated means of probabilities of effect being 1 given cause being either 1 for positive association or 0 for negative association using CondProb() and bootstrapping to estimate.}
#' \item{i}{An index}
#' \item{j}{An index}
#' \item{causalInfo\code{[['i,j']]}$CDirConfValInv}{An \code{alpha}*100th percentile confidence interval of estimated conditional probability of effect j being 1/0 given cause i's value being either the same (positive association) or opposite (negative association).}
#' \item{causalInfo\code{[['i,j']]}$CDirConfInv}{An \code{alpha}*100th percentile confidence interval of estimated causal direction degree of i cause j. }
#' \item{causalInfo\code{[['i,j']]}$CDirmean}{A mean-estimated-causal-direction degree of i cause j.}
#' \item{causalInfo\code{[['i,j']]}$testRes2}{A Mann-Whitney hypothesis test result for existence of causal direction. The null hypothesis is that the distributions of causal-direction degrees of i,j differ by a location shift of \code{CausalThs} and the alternative is  that distributions of causal-direction degrees of i,j is shifted greater than \code{CausalThs}.}
#' \item{causalInfo\code{[['i,j']]}$testRes1}{A Mann-Whitney hypothesis test result for existence of association by odd differences from \code{oddDiffFunc()}. The null hypothesis is that the distributions of absolute odd difference of i,j differ by a location shift of \code{IndpThs} and the alternative is  that distributions of absolute odd difference of i,j is shifted greater than \code{IndpThs}.}
#' \item{causalInfo\code{[['i,j']]}$sign}{A direction of i,j association: 1 for positive, 0 for negative, and -1 for no association.}
#' \item{causalInfo\code{[['i,j']]}$SignConfInv}{An \code{alpha}*100th percentile confidence interval of i,j odd difference from bootstrapping. }
#' \item{causalInfo\code{[['i,j']]}$Signmean}{A mean of i,j odd difference from bootstrapping.}
#'
#'@export
#'
bSCMCausalGraphFunc<-function(E1,Dboot,alpha=0.05,SignThs=0.05,CausalThs = 0.25)
{
  d<-dim(E1)[1]
  Ehat<-matrix(0,d,d)
  EValHat<-matrix(0,d,d)
  nboot= length(Dboot)
  bSignDist<-numeric(nboot)
  inxList<-c()

  causalInfo<-list()
  if(sum(E1) ==0)
    return(list(Ehat=E1,causalInfo=causalInfo))

  for(i in seq(1,d-1))
    for(j in seq(i+1,d))
    {
      if(E1[i,j]==1)
        inxList<-rbind(inxList, c(i,j))
    }

  signFlag<-numeric(dim(inxList)[1])-1
  # == Create the bootstrapping sequence of bSignDist

  for( itr in seq(dim(inxList)[1]) )
  {
    inx<-inxList[itr,]
    i=inx[1]
    j=inx[2]
    str<-sprintf("%d,%d",i,j)
    #print(str)
    for(k in seq(nboot))
    {
      bSignDist[k]<-oddDiffFunc(D=Dboot[[k]],i=inx[1],j=inx[2])
    }
    testRes1<-wilcox.test(x=abs(bSignDist), mu = SignThs, alternative = "greater")
    if(testRes1$p.value<=alpha)
    {
      if(mean(bSignDist) >0)
        signFlag[itr]<-1
      else
        signFlag[itr]<-0
    }

    if(signFlag[itr] != -1){
      #=========TODO dir inference
      bCausalDirDist<-numeric(nboot)

      bCausalDirValDistA<-numeric(nboot) # Y given Z
      bCausalDirValDistB<-numeric(nboot) # Z given Y
      for(k in seq(nboot))
      {
        D<-Dboot[[k]]

        z1<-numeric(d)-1
        y1<-numeric(d)-1
        y1[j]<-1
        z1[i]<-(signFlag[itr])
        a1<-CondProb(D,y=y1,z=z1)$condP
        b1<-CondProb(D,y=z1,z=y1)$condP
        bCausalDirDist[k]<-a1-b1
        bCausalDirValDistA[k]<-a1
        bCausalDirValDistB[k]<-b1
      }
      testRes2<-wilcox.test(x=abs(bCausalDirDist), mu = CausalThs, alternative = "greater")
      bmean <- mean( (bCausalDirDist) )
      #print(bmean)
      dirFlag=1
      if(testRes2$p.value<=alpha)
      {
        if(bmean >0) # i -> j
        {
          Ehat[i,j]<-1
          EValHat[i,j]<-mean(bCausalDirValDistA)
          causalInfo[[str]]$CDirConfValInv<-abs(quantile(bCausalDirValDistA, c(0+alpha/2, 1-alpha/2)) )
        }
        else #j -> i
        {
          Ehat[j,i]<-1
          EValHat[j,i]<-mean(bCausalDirValDistB)
          str<-sprintf("%d,%d",j,i)
          dirFlag=-1
          causalInfo[[str]]$CDirConfValInv<-abs(quantile(bCausalDirValDistB, c(0+alpha/2, 1-alpha/2)) )
        }
        causalInfo[[str]]$CDirConfInv<-abs(quantile(dirFlag*bCausalDirDist, c(0+alpha/2, 1-alpha/2)) )
        causalInfo[[str]]$CDirmean<-abs(bmean)
        causalInfo[[str]]$testRes2<-testRes2
        causalInfo[[str]]$testRes1<-testRes1
        causalInfo[[str]]$sign<-signFlag[itr]
        causalInfo[[str]]$SignConfInv<-quantile(bSignDist, c(0+alpha/2, 1-alpha/2))
        causalInfo[[str]]$Signmean<-mean(bSignDist)
      }
    }
  }
  return(list(Ehat=Ehat,causalInfo=causalInfo,EValHat=EValHat) )
}
