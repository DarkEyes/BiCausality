#'@title bSCMdeConfoundingGraphFunc function
#'
#'
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
            bIndpDist[b]<- max(c(indpFunc(nD,i,j,d,z1), indpFunc(nD,i,j,d,z2) ) )
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
