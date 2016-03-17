\name{mma-package}
\alias{mma-package}
\docType{package}
\title{
Mediation Analysis Package
}
\description{
This package is used to identify mediators and for general mediation analysis. Mediation effect refers to the effect conveyed by intervening variables to an observed relationship between an exposure and a response variable (outcome).  In this package, the exposure is called the predictor, the intervening variables are called mediators.  The mediation effects include the total effect, direct effect, and indirect effect, which are defined and the statistical inferences described in Yu et al.(2014).
}
\details{
\tabular{ll}{
Package: \tab mma\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-18\cr
License: \tab GPL (>= 2)\cr
}
\code{"\link[=data.org]{data.org}"} is used to identify potential mediators. \code{"\link[=med.binx]{med.binx}"} , \code{"\link[=boot.med.binx]{boot.med.binx}"}, \code{"\link[=med.contx]{med.contx}"}, and \code{"\link[=boot.med.contx]{boot.med.contx}"} are used for statistical inferences on the mediation effects when the predictor is binary or continuous according to the name of the function. \code{"\link[=mma]{mma}"} is a combined function that automatically identify potential mediators, based on which to make statistical inference on the mediation effects.
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}

Maintainer: Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi: 10.1037/0022-3514.51.6.1173>. "The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations," J. Pers Soc Psychol, 51(6), 1173-1182.

J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi: 10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi: 10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi: 10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\keyword{ Package }
\examples{
data("weight_behavior")
 #binary x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),mart=FALSE,nu=0.001,
   D=3,distn="bernoulli",family1=binomial(link = "logit"),n2=2)
 temp.b.b.mart<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),mart=TRUE,nu=0.05,
   D=3,family1=binomial(link = "logit"),n2=5)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 temp.b.c.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,jointm=list(n=1,j1=8:10),biny=FALSE, 
   predref="M",alpha=0.4,alpha2=0.4,n=2,seed=1,mart=FALSE,nu=0.05,D=3,
   distn="gaussian",family1=gaussian(link = "identity"),n2=2)
 temp.b.c.mart<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,jointm=list(n=1,j1=8:10, j2=12:13),biny=FALSE, 
   predref="M",alpha=0.4,alpha2=0.4, margin=1, n=2,seed=1,mart=TRUE,
   nu=0.05,distn="gaussian",family1=gaussian(link = "identity"),n2=2)
}
