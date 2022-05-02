BiCausality: Binary Causality Inference Framework
===========================================================
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-orange.svg)](https://spdx.org/licenses/MIT.html)

A framework to infer causality on binary data using frequent pattern minging and estimation statistics. Given a set of individual vectors S={x} where x(i) is a realization value of binary variable i, the framework infer empirical causal relations of binary variables i,j from S in a form of causal graph G=(V,E) where V is a set of nodes representing binary variables and there is an edge from i to j in E if the variable i causes j. The framework determines dependency among variables as well as analyzing confounding factors before decide whether i causes j. 

Note: The causal relations inferred by this work is not the real causal relations; they are empirical causal relations that needed
to be validated. Our main goal is to develop an exploratory data analysis tools to pinpoint possible causal relations to support
researchers before the validation in the field studies to find real causal relations. 

Installation
------------

For the newest version on github, please call the following command in R terminal.


``` r
remotes::install_github("DarkEyes/BiCausality")
```
This requires a user to install the "remotes" package before installing BiCausality.


Example: Inferred binary causal graph from simulation
----------------------------------------------------------------------------------
In the first step, we generate a simulation dataset as an input.
``` r
seedN<-2022

n<-200 # 200 individuals
d<-10 # 10 variables
mat<-matrix(nrow=n,ncol=d) # the input of framework

#Simulate binary data from binomial distribution where the probability of value being 1 is 0.5.
for(i in seq(n))
{
  set.seed(seedN+i)
  mat[i,] <- rbinom(n=d, size=1, prob=0.5)
}

mat[,1]<-mat[,2] | mat[,3]  # 1 causes by 2 and 3
mat[,4] <-mat[,2] | mat[,5] # 4 causses by 2 and 5
mat[,6] <- mat[,1] | mat[,4] # 6 causes by 1 and 4

```

We use the following function to infer whether X causes Y.
```{r}
# Run the function
library(BiCausality)
resC<-BiCausality::CausalGraphInferMainFunc(mat = mat,CausalThs=0.1, nboot =50, IndpThs=0.05)
```
The result of the ajacency matrix of the directed causal graph is below:

```r
resC$CausalGRes$Ehat
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 [1,]    0    0    0    0    0    1    0    0    0     0
 [2,]    1    0    0    1    0    0    0    0    0     0
 [3,]    1    0    0    0    0    0    0    0    0     0
 [4,]    0    0    0    0    0    1    0    0    0     0
 [5,]    0    0    0    1    0    0    0    0    0     0
 [6,]    0    0    0    0    0    0    0    0    0     0
 [7,]    0    0    0    0    0    0    0    0    0     0
 [8,]    0    0    0    0    0    0    0    0    0     0
 [9,]    0    0    0    0    0    0    0    0    0     0
[10,]    0    0    0    0    0    0    0    0    0     0
```
The value in the element EValHat[i,j] represents that i causes j if the value is not zero. For example, EValHat[2,1] = 1 implies node 2 causes node 1, which is correct since node 1 have nodes 2 and 3 as causal nodes.

The directed causal graph also can be plot using the code below.
```r
library(igraph)
net <- graph_from_adjacency_matrix(resC$CausalGRes$Ehat ,weighted = NULL)
plot(net, edge.arrow.size = 0.2, vertex.size =20 , vertex.color = '#D4C8E9',layout=layout_with_kk)
```
The plot is below.

<img src="https://github.com/DarkEyes/BiCausality/blob/master/man/FIG/causalGraph.png" width="550">


For the causal relation of variables 2 and 1, we can use the command below to see further information.

**Note that the odd difference between X and Y denoted oddDiff(X,Y) is define as
|P (X = 1, Y = 1) P (X = 0, Y = 0) −P (X = 0, Y = 1) P (X = 1, Y = 0)|.  If X is directly proportional to Y, then oddDiff(X,Y) is close to 1. If X is inverse of Y, then oddDiff(X,Y) is close to -1. If X and Y have no association, then oddDiff(X,Y) is close to zero.

```r
resC$CausalGRes$causalInfo[['2,1']]
```
Suppose Y is variable 1 and X is variable 2, the results are below.

```r
#This value represents the 95th percentile confidence interval of P(Y=1|X=1). 
$CDirConfValInv
 2.5% 97.5% 
    1     1 
#This value represents the 95th percentile confidence interval of |P(Y=1|X=1) - P(X=1|Y=1)|.
$CDirConfInv
     2.5%     97.5% 
0.3217322 0.4534494 

#This value represents the mean of |P(Y=1|X=1) - P(X=1|Y=1)|.
$CDirmean
[1] 0.3787904

#The test that has the null hypothesis that |P(Y=1|X=1) - P(X=1|Y=1)| below
#or equal the argument of parameter "CausalThs" and the alternative hypothesis
#is that |P(Y=1|X=1) - P(X=1|Y=1)| is greater than "CausalThs".
$testRes2

	Wilcoxon signed rank test with continuity correction

data:  abs(bCausalDirDist)
V = 1275, p-value = 3.893e-10
alternative hypothesis: true location is greater than 0.1


#The test that has the null hypothesis that |oddDiff(X,Y)| below 
#or equal the argument of parameter "IndpThs" and the alternative hypothesis is
#that |oddDiff(X,Y)| is greater than "IndpThs". 
$testRes1

	Wilcoxon signed rank test with continuity correction

data:  abs(bSignDist)
V = 1275, p-value = 3.894e-10
alternative hypothesis: true location is greater than 0.05

#If the test above rejects the null hypothesis with the significance threshold
#alpha (default alpha=0.05), then the value "sign=1", otherwise, it is zero.
$sign
[1] 1

#This value represents the 95th percentile confidence interval of oddDiff(X,Y)
$SignConfInv
      2.5%      97.5% 
0.08670325 0.13693900 

#This value represents the mean of oddDiff(X,Y)
$Signmean
[1] 0.1082242
```


Citation
----------------------------------------------------------------------------------
Chainarong Amornbunchornvej, Navaporn Surasvadi, Anon Plangprasopchok, and Suttipong Thajchayapong (2022). Framework for inferring empirical causal graphs
from binary variables to support multidimensional poverty analysis. **working on progress

Contact
----------------------------------------------------------------------------------
- Developer: C. Amornbunchornvej<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-3131-0370" href="https://orcid.org/0000-0003-3131-0370" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0003-3131-0370</a></div>
- <a href="https://www.nectec.or.th">Strategic Analytics Networks with Machine Learning and AI (SAI)</a>, <a href="https://www.nectec.or.th/en/">NECTEC</a>, Thailand
- Homepage: <a href="https://sites.google.com/view/amornbunchornvej/home">Link</a>
