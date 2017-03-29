\name{mma.par}
\alias{mma.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple Mediation Analysis Using Parallel Computing
}
\description{
Test for mediators and do statistical inferences on the identified mediation effects using parallel computing (multiple cores).
}
\usage{
mma.par(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,
                  catmed=NULL,catref=NULL,jointm=NULL,refy=NULL,predref=NULL,
                  alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
                  nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=NULL,n2=50,
                  weight=rep(1,nrow(x)),x.new=NULL,type=NULL,ncore=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates.
}
  \item{y}{
the vector of outcome variable.
}
  \item{pred}{
the column number that locates the predictor in x.
}
  \item{mediator}{
the list of mediators (column numbers in x or by variable names). The mediators to be checked can be identified by "contmed", "binmed" and "catmed", or by this argument, "mediator", where binary and categorical mediators in x are identified by factors, the reference group is the first level of the factor.  
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
  \item{refy}{
if y is binary, the reference group of y.By default, the reference group will be the first level of as.factor(y).
}
  \item{predref}{
if predictor is binary, identify the reference group of the binary predictor. By default, the reference group will be the first level of the predictor.
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
  \item{nonlinear}{
if TURE, Multiple Additive Regression Trees (MART) will be used to fit the final full model in estimating the outcome. The default value of nonlinear is FALSE, in which case, a generalized linear model will be used to fit the final full model.
}
  \item{df}{
if nonlinear is TURE, natural cubic spline will be used to fit the relationship between the countinuous predictor and each mediator.  The df is the degree of freedom in the ns() function, the default is 1.
}
  \item{nu}{
set the parameter "interaction.depth" in gbm function if MART is to be used, by default, nu=0.001. See also help(gbm.fit).
}
  \item{D}{
set the parameter "shrinkage" in gbm function if MART is to be used, by default, D=3. See also help(gbm.fit).
}
  \item{distn}{
the assumed distribution of the outcome if MART is used for final full model. The default value of distn is "bernoulli" for binary y, and "gaussian" for continuous y.
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is binomial(link = "logit") for binary y, gaussian(link="identity") for continuous y.
}
  \item{n2}{
the number of times of bootstrap resampling.  The default value is n2=50.
}
  \item{weight}{
the weight for each observation.
}
  \item{x.new}{
for continuous predictor only. of the same format as x (after data.org), with a new set of predictor and corresponding covariates on which to calculate the mediation effects. If is NULL, the mediation effects will be calculated based on the original data set.
}
  \item{type}{
the type of prediction when y is class Surv. Is "risk" if not specified.
}
  \item{ncore}{
the number of threads to be used for parallel computing.
}
}
\details{
mma.par is the same as mma, but it uses multiple cores for parallel computing to save the computing time.  
}
\value{
Returns an mma object.
  \item{estimation }{list the estimation of ie (indirect effect), te  (total effect), and de (direct effect from the predictor) separately.}  \item{bootsresults }{a list where the first item, ie, is a matrix of n2 rows where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) from the n2 bootstrap samples; the second item, te, is a vector of estimated total effects from the bootstrap sample; and the 3rd item, de, is a vector of estimated direct effect of the predictor from the bootstrap sample.} 
  \item{model }{a list where the first item, MART, is T if mart is fitted for the final model; the second item, Survival, is T if a survival model is fitted; the third item, type, is the type of prediction when a survival model is fitted; the fourth item, model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the fourth item, best.iter is the number of best iterations if MART is used to fit the final model.}
  \item{data }{a list that contains all the used data: x=x, y=y, dirx=dirx, binm=binm, contm=contm, catm=catm, jointm=jointm, binpred=F.}
  \item{boot.detail }{for continuous predictor only: a list that contains the mediation effects on each row of new.x: new.x=new.x, te1, de1, ie1.}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi:10.1037/0022-3514.51.6.1173>. "The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations," J. Pers Soc Psychol, 51(6), 1173-1182.

J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi:10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi:10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}.
}
\seealso{
\code{"\link[=data.org]{data.org}"} is for mediator tests, and \code{"\link[=med]{med}"} , and \code{"\link[=boot.med]{boot.med}"} for mediation analysis and inferences.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Mediator Tests }% __ONLY ONE__ keyword per line
