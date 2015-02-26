\name{boot.med.contx}
\alias{boot.med.contx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical Inference on Mediation Analysis with Continuous Predictor
}
\description{
To make inferences on the mediation effects when the predictor is continuous.
}
\usage{
boot.med.contx(x, y, dirx, binm = NULL, contm = NULL, catm = NULL, 
 jointm = NULL, margin = 1, n = 20, seed = sample(1:1000, 1), 
 mart = F, nu = 0.001, D = 3, distn = "gaussian", 
 family1 = gaussian(link = "identity"), n2 = 50)
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
  \item{binm}{
the column number of x that locates the binary mediators.
}
  \item{contm }{the column numbers of x that locate the potential continuous mediators.}
  \item{catm }{categorical mediators should be binarized and be presented as a list, where the first item is the number of categorical variables and the following items are the column numbers of each binarized categorical variable in x. data.org and data.org2 organize the categorical mediators in this format after they pass the mediator tests.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the column numbers of the mediators in x for each group of joint mediators.}
  \item{margin}{
the change in predictor when calculating the mediation effects, see Yu et al. (2014).
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
the assumed distribution of the outcome if MART is used for final full model. The default value of distn is "gaussian". If y is binary, may use distn="bernoulli".
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is gaussian(link="identity"). If y is binary, family1 can be defined as binomial(link = "logit").
}
  \item{n2}{
the number of times of bootstrap resampling.  The default value is n2=50.
}
}
\details{
The mediators are not tested in this function.  data.org or data.org2 should be used for the tests and the resulted list can be used directly to define the parameters in this function. boot.med.contx considers all variables in x as mediators or covariates in the final model and all variables identified by contm, binm, catm, or jointm as mediators.
}
\value{
Returns an mma object.
  \item{estimation }{list the estimation of ie (indirect effect), te  (total effect), and de (direct effect from the predictor) separately.}  \item{bootsresults }{a list where the first item, ie, is a matrix of n2 rows where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) from the n2 bootstrap samples; the second item, te, is a vector of estimated total effects from the bootstrap sample; and the 3rd item, de, is a vector of estimated direct effect of the predictor from the bootstrap sample.} 
  \item{model }{a list where the first item, mart, is T if MART is fitted for the final model; the second item, model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the third item, best.iter is the number of best iterations if MART is used to fit the final model.}
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014). "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}}

\seealso{
\code{"\link[=med.contx]{med.contx}"} to estimate mediation effects,  \code{"\link[=med.binx]{med.binx}"} and \code{"\link[=boot.med.binx]{boot.med.binx}"} for mediation analysis when the preditor is binary.
}
\examples{
data("weight_behavior")
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.contx<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,binpred=FALSE,alpha=0.4,alpha2=0.4)
 temp5<-boot.med.contx(x=data.contx$x,y,dirx=data.contx$dirx,
   binm=data.contx$binm, contm=data.contx$contm,catm=data.contx$catm,
   seed=1,n=2,family1=binomial(link = "logit"),n2=2)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.contx<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,binpred=FALSE,
   family1=gaussian(link="identity"),alpha=0.4,alpha2=0.4)
 temp7<-boot.med.contx(x=data.contx$x,y,dirx=data.contx$dirx,
   binm=data.contx$binm, contm=data.contx$contm,catm=data.contx$catm,
   seed=1,n=2,n2=2) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Continuous Predictor }% __ONLY ONE__ keyword per line
