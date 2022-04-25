#'@title bSCMCausalGraphFunc function
#'
#'
#'
#'@export
#'
bSCMCausalGraphFunc<-function(E1,Dboot,alpha=0.05,SignThs=0.05,CausalThs = 0.25,slack=0.001)
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
      bSignDist[k]<-oddDiffFunc(D=Dboot[[k]],i=inx[1],j=inx[2],slack=slack)
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
