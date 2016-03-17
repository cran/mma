\name{mma}
\alias{mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple Mediation Analysis
}
\description{
Test for mediators and do statistical inferences on the identified mediation effects.
}
\usage{
mma(x, y, pred, contmed = NULL, binmed = NULL, 
 binref = rep(1,length(binmed)), catmed = NULL, 
 catref = rep(1,length(catmed)), jointm = NULL, biny = T, binpred = T, 
 predref = 1, alpha = 0.1, alpha2 = 0.1, margin = 1, n = 20, 
 seed = sample(1:1000, 1), mart = F, nu = 0.001, D = 3,  
 distn = "bernoulli", family1 = binomial(link = "logit"), n2 = 50,
 weight=rep(1,length(y)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates
}
  \item{y}{
the vector of outcome variable.
}
  \item{pred}{
the column number that locates the predictor in x.
}
  \item{contmed}{
a vector of column numbers that locate the potential continuous mediators in x.
}
  \item{binmed}{
a vector of column numbers that locate the potential binary mediators in x.
}
  \item{binref}{
the defined reference groups of the binary potential mediators in binmed.
}
  \item{catmed}{
a vector of column numbers that locate the potential categorical mediators in x.
}
  \item{catref}{
the defined reference groups of the categorical potential mediators in catmed.
}
  \item{jointm}{
a list that identifies the mediators that need to be considered jointly, where the first item indicates the number of groups of mediators to be considered jointly, and each of the following items identifies the column numbers of the mediators in x for each group of joint mediators.
}
  \item{biny}{
if TRUE, the outcome variable is binary.  The default value of biny is TRUE.
}
  \item{binpred}{
if TRUE, the predictor is binary.  The default value of binpred is TRUE.
}
  \item{predref}{
if binpred is TRUE, identify the reference group of the binary predictor.
}
  \item{alpha}{
the significance level at which to test if the potential mediators (identified by contmed, binmed, and catmed) can be used as a covariate or mediator in estimating y when all variables in x are included in the model.  The default value is alpha=0.1
}
  \item{alpha2}{
the significant level at which to test if a potential mediator is related with the predictor. The default value is alpha2=0.1.
}  
  \item{margin}{
if binpred is FALSE, define the change in predictor when calculating the mediation effects, see Yu et al. (2014).
}
  \item{n}{
the time of resampling in calculating the indirect effects, default is n=20, see Yu et al. (2014).
}
  \item{seed}{
set seed to make the calculation repeatable. The default value of seed is from sample(1:1000,1).
}
  \item{mart}{
if TURE, Multiple Additive Regression Trees (MART) will be used to fit the final full model in estimating the outcome. The default value of mart is FALSE, in which case, a generalized linear model will be used to fit the final full model.
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
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is binomial(link = "logit"). If y is binary, family1 can be defined as gaussian(link="identity").
}
  \item{n2}{
the number of times of bootstrap resampling.  The default value is n2=50.
}
  \item{weight}{
the weight for each observation.
}
}
\details{
mma first tests if the potential mediators defined by binm, contm, and catm should be treated as mediators or covariates (if none, the variable will be deleted from further analysis). All variables identified by jointm are treated as mediators. All other variables in x that are not tested are treated as covariates. Then mma does the mediation effects estimation and inference on the selected variables.  
}
\value{
Returns an mma object.
  \item{estimation }{list the estimation of ie (indirect effect), te  (total effect), and de (direct effect from the predictor) separately.}  \item{bootsresults }{a list where the first item, ie, is a matrix of n2 rows where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) from the n2 bootstrap samples; the second item, te, is a vector of estimated total effects from the bootstrap sample; and the 3rd item, de, is a vector of estimated direct effect of the predictor from the bootstrap sample.} 
  \item{model }{a list where the first item, mart, is T if MART is fitted for the final model; the second item, model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the third item, best.iter is the number of best iterations if MART is used to fit the final model.}
  \item{data }{a list that contain all the used data}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi: 10.1037/0022-3514.51.6.1173>. "The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations," J. Pers Soc Psychol, 51(6), 1173-1182.

J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi: 10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi: 10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi: 10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\seealso{
\code{"\link[=data.org]{data.org}"} is for mediator tests, and \code{"\link[=med.binx]{med.binx}"} , \code{"\link[=boot.med.binx]{boot.med.binx}"}, \code{"\link[=med.contx]{med.contx}"}, and \code{"\link[=boot.med.contx]{boot.med.contx}"} for mediation analysis and inferences.
}
\examples{
data("weight_behavior")
 #binary x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),mart=FALSE,nu=0.001,
   D=3,n2=2)
 temp.b.b.mart<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),mart=TRUE,nu=0.05,
   D=3,n2=5)
 
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Mediator Tests }% __ONLY ONE__ keyword per line
