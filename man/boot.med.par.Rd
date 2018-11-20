\name{boot.med.par}
\alias{boot.med.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical Inference on Mediation Analysis Using Parallel Computing
}
\description{
To make inferences on the mediation effects when the predictor is continuous or binary using parallel computing.
}
\usage{
boot.med.par(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,
             contm=data$contm,catm=data$catm,jointm=data$jointm,cova=data$cova,margin=1,
             n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,
             D=3,distn=NULL,family1=data$family1,n2=50,w=rep(1,nrow(x)),
             refy=NULL,x.new=x,pred.new=dirx,binpred=data$binpred,
             type=NULL,w.new=NULL,ncore=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
the list of result from data.org that organize the covariates, mediators, predictor and outcome. If data is FALSE, then need to set x, y, dirx, contm, catm, and jointm.
}
   \item{x}{
a data frame contains all mediators and covariates.  Need to set up only when data is FALSE. All varaibles in x that are not identified as potential mediators (by mediator, contmed, binmed, catmed, or jointm) are forced in mediation analysis as covariates.
}
  \item{y}{
the vector of outcome variable.  Need to set up only when data is FALSE.
}
  \item{dirx }{the vector/matrix of predictor(s). Need to set up only when data is FALSE.}
  \item{binm}{
the variable names or the column number of x that locates the binary mediators. Need to set up only when data is FALSE.
}
  \item{contm }{the variable names or the column numbers of x that locate the potential continuous mediators. Need to set up only when data is FALSE.}
  \item{catm }{categorical mediators should be binarized and be presented as a list, where the first item is the number of categorical variables and the following items are the names or the column numbers of each binarized categorical variable in x. data.org organizes the categorical mediators in this format after they pass the mediator tests. Need to set up only when data is FALSE.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the names or the column numbers of the mediators in x for each group of joint mediators. Need to set up only when data is FALSE.}
  \item{cova}{the covariate names in x for predicting the mediators. The cova cannot be potential mediators, so cova are the covariates for the response automatically. If the covariates are for all potential mediators, cova is a vector of covaraite names. If the covariates are for some specific potential mediators, cova is a list with the first item be the names of potential mediators that use the covariates, and the second item be the names of covariates. Works mainly for continuous predictors (pred), also the specified covariates should have no missing if only some mediators uses the covariates.
}
  \item{margin}{
the change in predictor when calculating the mediation effects, see Yu et al. (2014).
}
  \item{n}{
the time of resampling in calculating the indirect effects, default is n=20, see Yu et al. (2014).
}
  \item{seed}{
set seed to make the calculation repeatable. The default value of seed is from sample(1:1000,1).
}
  \item{nonlinear}{
if TURE, Multiple Additive Regression Trees (MART) will be used to fit the final full model in estimating the outcome. The default value of nonlinear is FALSE, in which case, a generalized linear model will be used to fit the final full model.
}
  \item{df}{
if nonlinear is TURE, natural cubic spline will be used to fit the relationship between the predictor and each mediator.  The df is the degree of freedom in the ns() function, the default is 1.
}
  \item{nu}{
set the parameter "interaction.depth" in gbm function if MART is to be used, by default, nu=0.001. See also help(gbm.fit).
}
  \item{D}{
set the parameter "shrinkage" in gbm function if MART is to be used, by default, D=3. See also help(gbm.fit).
}
  \item{distn}{
the assumed distribution of the outcome if MART is used for final full model. The default value of distn is "gaussian" for coninuous y, and "bernoulli" for binary y.
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is gaussian(link="identity") for continuous y, and binomial(link = "logit") for binary y.
}
  \item{n2}{
the number of times of bootstrap resampling.  The default value is n2=50.
}
  \item{w}{
the weight for observations.
}
  \item{refy}{
if y is binary, the reference group of y.
}
  \item{binpred}{
if TRUE, the predictor is binary.
}
  \item{type}{
the type of prediction when y is class Surv. Is "risk" if not specified.
}
  \item{x.new}{
of the same format as x, with a new set of covariates and mediators on which to calculate the mediation effects.
}
  \item{pred.new}{
a new set of predictor(s).
}
  \item{w.new}{
the weights for new.x.
}
  \item{ncore}{
the number of threads to be used for parallel computing.
}
}
\details{
Same as boot.med, except that this function uses parallel computing. Need to call the package "doParallel" for this function. 
}
\value{
Returns an mma object.
  \item{estimation }{list the estimation of ie (indirect effect), te  (total effect), and de (direct effect from the predictor) separately.}  \item{bootsresults }{a list where the first item, ie, is a matrix of n2 rows where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) from the n2 bootstrap samples; the second item, te, is a vector of estimated total effects from the bootstrap sample; and the 3rd item, de, is a vector of estimated direct effect of the predictor from the bootstrap sample.} 
  \item{model }{a list where the first item, MART, is T if mart is fitted for the final model; the second item, Survival, is T if a survival model is fitted; the third item, type, is the type of prediction when a survival model is fitted; the fourth item, model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the fourth item, best.iter is the number of best iterations if MART is used to fit the final model.}
  \item{data }{a list that contains all the used data: x=x, y=y, dirx=dirx, binm=binm, contm=contm, catm=catm, jointm=jointm, binpred=F.}
  \item{boot.detail }{a list that contains the mediation effects on each row of new.x: new.x=new.x, te1, de1, ie1.}
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}.}

\seealso{
\code{"\link[=med]{med}"} just estimate the mediation effects.
}
\examples{
\donttest{
library(doParallel)
data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-boot.med.par(data=data.bin,n=2,n2=4,seed=1)
 temp2<-boot.med.par(data=data.bin,n=2,n2=4,nu=0.05,seed=1,nonlinear=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Continuous Predictor }% __ONLY ONE__ keyword per line
