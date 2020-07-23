#'@title CausalGraphInferMainFunc function
#'
#'
#'
#'@export
#'
CausalGraphInferMainFunc<-function(mat,alpha=0.05,nboot=100,IndpThs=0.05,CausalThs = 0.25,slack=0.001)
{
  print("Inferring dependent graph")
  res<-bSCMDepndentGraphFunc(mat,nboot=nboot,alpha=alpha,IndpThs=IndpThs)
  print("Removing confounder(s)")
  res2<-bSCMdeConfoundingGraphFunc(res,IndpThs=IndpThs,alpha=alpha)
  print("Inferring causal graph")
  res3<-bSCMCausalGraphFunc(res2$E1,res$Dboot ,alpha=alpha,SignThs=IndpThs,CausalThs = CausalThs,slack=slack)
  return(list(depRes = res, ConfoundRes= res2, CausalGRes= res3))
}
