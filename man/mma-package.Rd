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
Version: \tab 6.1\cr
Date: \tab 2018-03-08\cr
License: \tab GPL (>= 2)\cr
}
\code{"\link[=data.org]{data.org}"} is used to identify potential mediators. \code{"\link[=med]{med}"} , and \code{"\link[=boot.med]{boot.med}"} are used for statistical inferences on the mediation effects when the predictor is binary or continuous. \code{"\link[=mma]{mma}"} is a combined function that automatically identify potential mediators, based on which to make statistical inference on the mediation effects.
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}

Maintainer: Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi:10.1037/0022-3514.51.6.1173>. "The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations," J. Pers Soc Psychol, 51(6), 1173-1182.

J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi:10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi:10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Yu, Q., Wu, X., Li, B., and Scribner, R. (2018). <doi:10.1002/sim.7977>. "Multiple Mediation Analysis with Survival Outcomes – With an Application to Explore Racial Disparity in Breast Cancer Survival," Statistics in Medicine.

Yu, Q., Medeiros, KL, Wu, X., and Jensen, R. (2018). <doi:10.1007/s11336-018-9612-2>. "Explore Ethnic Disparities in Anxiety and Depression Among Cancer Survivors Using Nonlinear Mediation Analysis," Psychometrika, 83(4), 991-1006. 
}
\keyword{ Package }
\examples{
data("weight_behavior")
#binary predictor
 #binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 \donttest{
 temp.b.b.mart<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,nonlinear=TRUE,n=2,n2=5)
 #continuous y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,1]
 temp.b.c.glm<-mma(x,y,pred=pred,mediator=5:12,jointm=list(n=1,j1=7:9), 
                     predref="M",alpha=0.4,alpha2=0.4,n2=20)
 temp.b.c.mart<-mma(x,y,pred=pred,mediator=5:12,jointm=list(n=1,j1=7:9), 
                     predref="M",alpha=0.4,alpha2=0.4,
                     n=2,seed=1,nonlinear=TRUE,n2=20)}
}
