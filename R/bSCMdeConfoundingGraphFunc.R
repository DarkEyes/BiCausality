#'@title bSCMdeConfoundingGraphFunc function
#'
#' @description
#' This function removes any association/dependency of variables i,j that have any confounding factor k s.t. given k, i and j are independent.
#'
#' @param dat is the result of inferring dependencies between all pairs of variables from \code{bSCMDepndentGraphFunc()}.
#' @param alpha is a significance threshold for hypothesis tests (Mann Whitney)
#'  that deploys for testing degrees of dependency, association direction, and causal direction. The default is 0.5.
#' @param IndpThs is a threshold for the degree of dependency. In the independence test, to claim that any variables are dependent, the dependency degree must greater than this value significantly. The default is 0.05.
#'
#' @return This function returns an adjacency matrix of dependencies that have no confounding factors.
#' \item{E1}{An adjCausalGRes$causalInfo[['i,j']]acency matrix of undirected graph after filtering associations without true causal directions from any confounding factor.}
#' \item{E2}{A matrix of associations that have confounding factors where \code{E2[i,j]=0} if no confounding factor and \code{E2[i,j]=k} if k is a confounding factor of i and j.}
#'
#' @examples
#' #bSCMdeConfoundingGraphFunc(simData$resC$depRes)
#'
#'@export
#'
bSCMdeConfoundingGraphFunc<-function(dat,IndpThs=0.05,alpha=0.05)
{
  E0<-dat$E0
  d<-dim(E0)[1]
  E1<-matrix(0,d,d)
  E2<-matrix(0,d,d)
  nboot<-length(dat$Dboot)
  for(info in dat$depInfo)
  {
    i<-info$indices[1]
    j<-info$indices[2]
    # i and j have another dependent var
    if(E0[i,j]==1 && sum(E0[i,])>1 && sum(E0[j,])>1)
    {
      z<- 1:d
      z<-z[E0[i,] & E0[j,]]
      flag=0
      if(length(z)>0)
        for(z0 in z)
        {
          #print(sprintf("%d,%d | %d",i,j,z0))
          bIndpDist<-numeric(nboot)
          z1<-numeric(d)-1
          z1[z0]<-1
          z2<-numeric(d)-1
          z2[z0]<-0
          for(b in seq(nboot) )
          {
            nD<-dat$Dboot[[b]]
            bIndpDist[b]<- max(c(indpFunc(nD,i,j,z1), indpFunc(nD,i,j,z2) ) )
          }
          testRes<-wilcox.test(x=bIndpDist, mu = IndpThs, alternative = "greater")
          #print(sprintf("pval:%f",testRes$p.value))
          #check whether i is dependent with j given z0
          if(testRes$p.value>alpha)
          {
            flag=1
            E2[i,j]<-z0
            break
          }
        }
      if(flag==0)
      {
        E1[i,j]<-1
        E1[j,i]<-1
      }

    }
    else
    {
      E1[i,j]<-1
      E1[j,i]<-1
    }
  }
  return(list(E1=E1,E2=E2))
}
