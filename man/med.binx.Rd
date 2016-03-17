\name{med.binx}
\alias{med.binx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mediation Analysis with Binary Predictor
}
\description{
To estimate the mediation effects when the predictor is binary. 
}
\usage{
med.binx(x, y, dirx, contm = NULL, catm = NULL, jointm = NULL, 
  allm = c(contm, catm), n = 20, seed = sample(1:1000, 1), mart = F,
  nu = 0.001, D = 3, distn = "bernoulli", family1 = binomial("logit"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{
a data frame contains the predictor, all mediators and covariates.
}
  \item{y}{
the vector of outcome variable.
}
  \item{dirx }{the column number of x that locates the predictor.}
  \item{contm }{the column numbers of x that locate the potential continuous mediators.}
  \item{catm }{the column numbers of x that locate the potential binary or categorical mediators.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the column numbers of the mediators in x for each group of joint mediators.}
  \item{allm}{
the column numbers of all mediators.  The default value of allm is c(contm,catm).
}
  \item{n}{
the time of resampling when calculating the indirect effects, see Yu et al. (2014). The default value is n=20.
}
  \item{seed}{
set seed to make the calculation repeatable. The default value of seed is sampled by sample(1:1000,1).
}
  \item{mart}{
if TURE, Multiple Additive Regression Trees (MART) will be used to fit the final full model in estimating the outcome. The default value of mart is FALSE, in which case a generalized linear model will be used to fit the final full model.
}
  \item{nu}{
set the parameter "interaction.depth" in gbm function if MART is to be used, by default, nu=0.001. See also help(gbm.fit).
}
  \item{D}{
 set the parameter "shrinkage" in gbm function if MART is to be used, by default, D=3. See also help(gbm.fit).
}
  \item{distn}{
the assumed distribution of the outcome if MART is used for final full model. The default value of distn is "bernoulli". If y is continuous, may use distn="gaussian".
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is binomial(link = "logit"). If y is continuous, family1 can be defined as gaussian(link="identity").
}
}
\details{
The mediators are not tested in this function.  data.org should be used for the tests and the resulted list can be used directly to define the parameters in this function. med.binx considers all variables in x as mediators or covariates in the final model and all variables identified by contm, catm, allm, or jointm as mediators.
}
\value{
  \item{denm }{a matrix where each column gives the estimated direct effect not from the corresponding mediator (identified by the column name) and each row corresponding to the results from one resampling.}
  \item{ie }{a matrix where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) and each row corresponding to the results from one resampling.} 
  \item{te }{a vector of the estimated total effect from each resampling.}
  \item{model }{the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables.}
  \item{best.iter }{the number of best iterations if MART is used to fit the final model, is NULL if the final model is a generalized linear model.}
}
\references{
J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi: 10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi: 10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi: 10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}
\seealso{
\code{"\link[=boot.med.binx]{boot.med.binx}"} to make inferences on the estimated mediation effects using bootstrap method, \code{"\link[=med.contx]{med.contx}"}, and \code{"\link[=boot.med.contx]{boot.med.contx}"} for mediation analysis when the preditor is continuous.
}
\examples{
data("weight_behavior")
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-med.binx(x=data.bin$x,y,dirx=data.bin$dirx,
   contm=data.bin$contm,catm=data.bin$catm,n=2,seed=1)
 temp2<-med.binx(x=data.bin$x,y,dirx=data.bin$dirx,
   contm=data.bin$contm,catm=data.bin$catm,n=2,seed=1,mart=TRUE)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.cont<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,biny=FALSE, 
   family1=gaussian(link="identity"),predref="M",alpha=0.4,alpha2=0.4)
 temp3<-med.binx(data.cont$x,y,dirx=data.cont$dirx,contm=data.cont$contm,
   catm=data.cont$catm,n=2,seed=1,family1=gaussian(link="identity")) 
 temp4<-med.binx(data.cont$x,y,dirx=data.cont$dirx,contm=data.cont$contm,
   catm=data.cont$catm,n=2,seed=1,mart=TRUE,distn="gaussian") 
}
\keyword{ Mediation Analysis }
\keyword{ Binary Predictor }
