\name{med.par}
\alias{med.par}
\title{
Mediation Analysis Using Parallel Computing
}
\description{
To estimate the mediation effects when the predictor is binary or continuous.
}
\usage{
med.par(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm,
        contm = data$contm, catm = data$catm, jointm = data$jointm, 
        allm = c(contm, catm), margin=1,n=20,seed=sample(1:1000,1), 
        nonlinear=F, df=1, nu=0.001,D=3,distn=NULL,family1=NULL,refy=0,
        binpred=data$binpred,x.new=x,type=NULL,ncore=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
the list of result from data.org that organize the covariates, mediators, predictor and outcome. If data is FALSE, then need to set x1, y1, dirx, contm, catm, and jointm.
}
   \item{x}{
a data frame contains the predictor, all mediators and covariates.  Need to set up only when data is FALSE.
}
  \item{y}{
the vector of outcome variable.  Need to set up only when data is FALSE.
}
  \item{dirx }{the predictor name or the column number of x that locates the predictor. Need to set up only when data is FALSE.}
  \item{binm}{
the variable names or the column number of x that locates the binary mediators. Need to set up only when data is FALSE.
}
  \item{contm }{the variable names or the column numbers of x that locate the potential continuous mediators. Need to set up only when data is FALSE.}
  \item{catm }{categorical mediators should be binarized and be presented as a list, where the first item is the number of categorical variables and the following items are the names or the column numbers of each binarized categorical variable in x. data.org organizes the categorical mediators in this format after they pass the mediator tests. Need to set up only when data is FALSE.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the names or the column numbers of the mediators in x for each group of joint mediators. Need to set up only when data is FALSE.}
  \item{allm}{
the column numbers of all mediators.  Need to set up only when data is FALSE. The default value of allm is c(contm,catm).
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
the assumed distribution of the outcome if MART is used for final full model. The default value of distn is "gaussian" for continuous y and "bernoulli" for binary y.
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is gaussian(link="identity") for contiuous y, and binomial(link = "logit") for binary y.
}
  \item{refy}{
if y is binary, the reference group of y.
}
  \item{binpred}{
if TRUE, the predict variable is binary.  
}
  \item{x.new}{
of the same format as x, with a new set of predictor and corresponding covariates on which to calculate the mediation effects. For continuous predictor only.
}
  \item{type}{
the type of prediction when y is class Surv. By default, type is "risk".
}
  \item{ncore}{
the number of threads to be used for parallel computing.
}
}
\details{
The med.par function is the med function using paralle computing.
}
\value{The result is an med object with:
  \item{denm }{a matrix where each column gives the estimated direct effect not from the corresponding mediator (identified by the column name), see Yu et al. (2014) for the definition, and each row corresponding to the results from one resampling for binary predictor or the results on a row of x.new for continuous predictor.}
  \item{ie }{a matrix where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) and each row corresponding to the results from one resampling for binary predictor or the results on a row of x.new for continuous predictor.} 
  \item{te }{a vector of the estimated total effect on x.new.}
  \item{model }{a list, where the first item, MART, is TRUE if a mart is fitted as the final model; the second item, Survial is T if a survival model is fit; the third item, type, is the type of prediction when a survival model is fitted; the fourth item, full.model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the fifth item, best.iter,  is the number of best iterations if MART is used to fit the final model, is NULL if the final model is a generalized linear model.}
}
\references{
J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi:10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi:10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}.
}

\seealso{
\code{"\link[=boot.med]{boot.med}"} to make inferences on the estimated mediation effects using bootstrap method.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
