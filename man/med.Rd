\name{med}
\alias{med}
\title{
Mediation Analysis with Binary or Continuous Predictor
}
\description{
To estimate the mediation effects when the predictor is binary or continuous.
}
\usage{
med(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm, 
    contm = data$contm,catm = data$catm, jointm = data$jointm, cova=data$cova,
    allm = c(contm, catm), margin=1,n=20,seed=sample(1:1000,1), 
    nonlinear=F, df=1, nu=0.001,D=3,distn=NULL,family1=data$family1,
    refy=rep(0,ncol(y)),binpred=data$binpred,x.new=x,pred.new=dirx,
    type=NULL,w=NULL,w.new=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
the list of result from data.org that organize the covariates, mediators, predictor and outcome. If data is FALSE, then need to set x1, y1, dirx, contm, catm, and jointm.
}
   \item{x}{
a data frame contains all mediators and covariates.  Need to set up only when data is FALSE. All varaibles in x that are not identified as potential mediators (by mediator, contmed, binmed, catmed, or jointm) are forced in mediation analysis as covariates.
}
  \item{y}{
the vector of outcome variable.  Need to set up only when data is FALSE.
}
  \item{dirx }{the vector or matrix of predictor(s). The reference group is set to be 0. Need to set up only when data is FALSE.}
  \item{binm}{
the variable names or the column number of x that locates the binary mediators. Need to set up only when data is FALSE.
}
  \item{contm }{the variable names or the column numbers of x that locate the potential continuous mediators. Need to set up only when data is FALSE.}
  \item{catm }{categorical mediators should be binarized and be presented as a list, where the first item is the number of categorical variables and the following items are the names or the column numbers of each binarized categorical variable in x. data.org organizes the categorical mediators in this format after they pass the mediator tests. Need to set up only when data is FALSE.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the names or the column numbers of the mediators in x for each group of joint mediators. Need to set up only when data is FALSE.}
  \item{cova}{the covariate names in x for predicting the mediators. The cova cannot be potential mediators, so cova are the covariates for the response automatically. If the covariates are for all potential mediators, cova is a vector of covaraite names. If the covariates are for some specific potential mediators, cova is a list with the first item be the names of potential mediators that use the covariates, and the second item be the names of covariates. Need to set up only when data is FALSE.
}
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
the assumed distribution(s) of the outcome(s) if MART is used for final full model. The default value of distn is "gaussian" for continuous y,  "bernoulli" for binary y and coxph for "Surv" class outcome.
}
  \item{family1}{a list with the ith item 
define the conditional distribution of y[,i] given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is gaussian(link="identity") for contiuous y[,i], and binomial(link = "logit") for binary y[,i].
}
  \item{refy}{
if y is binary, the reference group of y.
}
  \item{binpred}{
if TRUE, the predict variable is binary.  
}
  \item{x.new}{
A new set of predictor and corresponding covariates, of the same format as x (after data.org),  on which to calculate the mediation effects. For continuous predictor only. If is NULL, the mediation effects will be calculated based on the original data set.
}
  \item{pred.new}{
A new set of predictor(s), of the same format as x (after data.org), on which to calculate the mediation effects. For continuous predictor only.
}
  \item{type}{
the type of prediction when y is class Surv. By default, type is "risk".
}
  \item{w}{
the weight for each case in x.
}
  \item{w.new}{
the weight for each case in x.new.
}
}
\details{
The mediators are not tested in this function.  data.org should be used first for the tests and data organizing, and then the resulting list from data.org can be used directly to define the arguments in this function. med considers all variables in x as mediators or covariates in the final model and all variables identified by contm, binm, catm, or jointm as mediators.
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
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

\seealso{
\code{"\link[=boot.med]{boot.med}"} to make inferences on the estimated mediation effects using bootstrap method.
}
\examples{
data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
  binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2,seed=1)
temp2<-med(data=data.bin,n=2,seed=1,nonlinear=TRUE)

\donttest{
 #multivariate predictor
 x=weight_behavior[,c(2:3,5:14)]
 pred=weight_behavior[,4]
 y=weight_behavior[,15]
 data.b.b.2.3<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 temp1.2<-med(data.b.b.2.3,n=2,seed=1)
 temp2.2<-med(data.b.b.2.3,n=2,seed=1,nonlinear=TRUE)

 #multivariate responses
 x=weight_behavior[,c(2:3,5:14)]
 pred=weight_behavior[,4]
 y=weight_behavior[,c(1,15)]
 data.b.b.2.4<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 temp1.3<-med(data.b.b.2.4,n=2,seed=1)
 temp2.3<-med(data.b.b.2.4,n=2,seed=1,nonlinear=TRUE)
}

#continuous y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,1]
 data.cont<-data.org(x,y,pred=pred,mediator=5:12,jointm=list(n=1,j1=7:9), 
                     predref="M",alpha=0.4,alpha2=0.4)
 temp3<-med(data=data.cont,n=2,seed=1) 
 temp4<-med(data=data.cont,n=2,seed=1,nonlinear=TRUE) 

##continuous x
#binary y
 x=weight_behavior[,3:14]
 pred=weight_behavior[,2]
 y=weight_behavior[,15]
 data.contx<-data.org(x,y,pred=pred,mediator=4:10,alpha=0.4,alpha2=0.4)
 temp5<-med(data=data.contx,seed=1,n=2)
 temp6<-med(data=data.contx,seed=1,n=2,nonlinear=TRUE,nu=0.05)

#continuous y
x=weight_behavior[,3:14]
y=weight_behavior[,1]
pred=weight_behavior[,2]
data.contx<-data.org(x,y,pred=pred,contmed=c(11:12),binmed=c(6,10),
                     binref=c(1,1),catmed=5,catref=1,
                     alpha=0.4,alpha2=0.4)
\donttest{
temp7<-med(data=data.contx,seed=1,n=2) 
temp8<-med(data=data.contx,seed=1,n=2,nonlinear=TRUE,nu=0.05) 
}
##Surv class outcome (survival analysis)
#data(cgd0)       #a dataset in the survival package
cgd1<-cgd0
x=cgd1[,c(4:5,7:12)]
pred=cgd1[,6]
status<-ifelse(is.na(cgd1$etime1),0,1)
y=Surv(cgd1$futime,status)          
#for continuous predictor
data.surv.contx<-data.org(x,y,pred=pred,mediator=1:ncol(x),      
                          alpha=0.5,alpha2=0.5)
temp9.contx<-med(data=data.surv.contx,seed=1,n=2,type="lp")

#close to mart results when use type="lp"
temp9.contx
temp10.contx<-med(data=data.surv.contx,n=2,seed=1,nonlinear=TRUE)  
#results in the linear part unit
temp10.contx

#for binary predictor
\donttest{
x=cgd1[,c(5:12)]
pred=cgd1[,4]
data.surv.binx<-data.org(x,y,pred=pred,mediator=1:ncol(x),   
                    alpha=0.4,alpha2=0.4)
temp9.binx<-med(data=data.surv.binx,seed=1,n=2,type="lp") 
temp9.binx
temp10.binx<-med(data=data.surv.binx,n=2,seed=1,nonlinear=TRUE)  
temp10.binx}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
