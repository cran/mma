\name{med}
\alias{med}
\title{
Mediation Analysis with Binary or Continuous Predictor
}
\description{
To estimate the mediation effects when the predictor is binary or continuous.
}
\usage{
med(data, x=data$bin.results$x, y=data$bin.results$y, dirx=data$bin.results$dirx, 
              binm=data$bin.results$binm,contm = data$bin.results$contm, 
              catm = data$bin.results$catm, jointm = data$bin.results$jointm, 
              cova=data$bin.results$cova, allm = c(contm, catm), 
              margin=1, n=20, nonlinear=F, df1=1, nu=0.001,D=3,distn=NULL,
              family1=data$bin.results$family1,refy=rep(0,ncol(y)),
              binpred=data$bin.results$binpred,x.new=x,pred.new=dirx, 
              cova.new=cova,type=NULL, w=NULL, w.new=NULL,xmod=NULL,custom.function=NULL)
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
  \item{cova}{The data frame of covariates to be used to predict the mediator(s). The covariates in cova cannot be potential mediators. If the covariates are for some specific potential mediators, cova$for.m is the vector of names of potential mediators that use the covariates. Works for continuous predictors (pred) only, also the specified covariates should have no missing if only some mediators uses the covariates.  Need to set up only when data is FALSE.
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
  %\item{seed}{
%set seed to make the calculation repeatable. The default value of seed is from sample(1:1000,1).
%}
  \item{nonlinear}{
if TURE, Multiple Additive Regression Trees (MART) will be used to fit the final full model in estimating the outcome. The default value of nonlinear is FALSE, in which case, a generalized linear model will be used to fit the final full model.
}
  \item{df1}{
if nonlinear is TURE, natural cubic spline will be used to fit the relationship between the predictor and each mediator.  The df is the degree of freedom in the ns() function, the default is 1.
}
  \item{nu}{
set the parameter "shrinkage" in gbm function if MART is to be used, by default, nu=0.001. See also help(gbm.fit).
}
  \item{D}{
set the parameter "interaction.depth" in gbm function if MART is to be used, by default, D=3. See also help(gbm.fit).
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
A new set of covariates, of the same format as x (after data.org),  on which to calculate the mediation effects. For continuous predictor only. If is NULL, the mediation effects will be calculated based on the original data set.
}
  \item{pred.new}{
A new set of predictor(s), of the same format as pred (after data.org), on which to calculate the mediation effects. For continuous predictor only.
}
  \item{cova.new}{
A new set of covariate(s) to estimate mediators, of the same format as cova.
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
  \item{xmod}{
If there is a moderator, xmod gives the moderator's name in cova and/or x.
}
  \item{custom.function}{a string of customer defined final model for predicting the outcome(s). The response variable should be noted as "responseY", the dataset should be noted as "dataset123". The weights for observations should be noted as "weights123". The covariates should be in x or pred. Use "~." for all varaibles in x and pred is allowed. The customer defined model should be able to make prediction by using "predict(object, newdata=...)", where the object is the results of the fitted model. If a specific package needs to be called to fit the model, the user should call the package first. For example, if the gamlss package is used to fit a piecewise polynomial with "age" to be the mediator and "race" be the predictor, we can set custom.function="gamlss(responseV~b(age,df=3)+race,data=dataset123,tace=FALSE)". The length of custom.function should be the same as the dimension of y. If custom.function[j] is NA, the usual method will be used to fit y[j]. 
}
}
\details{
The mediators are not tested in this function.  data.org should be used first for the tests and data organizing, and then the resulting list from data.org can be used directly to define the arguments in this function. med considers all variables in x as mediators or covariates in the final model and all variables identified by contm, binm, catm, or jointm as mediators.
}
\value{The result is an med object with the item a.binx for results of binary or categorical exposure(s) and the item a.contx for continuous exposure(s):
  \item{denm }{a matrix where each column gives the estimated direct effect not from the corresponding mediator (identified by the column name), see Yu et al. (2014) for the definition, and each row corresponding to the results from one resampling for binary predictor or the results on a row of x.new for continuous predictor.}
  \item{ie }{a matrix where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) and each row corresponding to the results from one resampling for binary predictor or the results on a row of x.new for continuous predictor.} 
  \item{te }{a vector of the estimated total effect on x.new.}
  \item{model }{a list, where the first item, MART, is TRUE if a mart is fitted as the final model; the second item, Survial is T if a survival model is fit; the third item, type, is the type of prediction when a survival model is fitted; the fourth item, full.model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the fifth item, best.iter,  is the number of best iterations if MART is used to fit the final model, is NULL if the final model is a generalized linear model.}
}
\references{
J.H. Friedman, T. Hastie, R. Tibshirani (2000) <doi:10.1214/aos/1016120463>. "Additive Logistic Regression: a Statistical View of Boosting," Annals of Statistics 28(2):337-374.

J.H. Friedman (2001) <doi:10.1214/aos/1013203451>. "Greedy Function Approximation: A Gradient Boosting Machine," Annals of Statistics 29(5):1189-1232.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Yu, Q., Wu, X., Li, B., and Scribner, R. (2018). <doi:10.1002/sim.7977>. "Multiple Mediation Analysis with Survival Outcomes – With an Application to Explore Racial Disparity in Breast Cancer Survival," Statistics in Medicine.

Yu, Q., Medeiros, KL, Wu, X., and Jensen, R. (2018). <doi:10.1007/s11336-018-9612-2>. "Explore Ethnic Disparities in Anxiety and Depression Among Cancer Survivors Using Nonlinear Mediation Analysis," Psychometrika, 83(4), 991-1006. 
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
temp1<-med(data=data.bin,n=2)
#or use self-defined final function
temp1<-med(data=data.bin,n=2,custom.function = 
           'glm(responseY~.,data=dataset123,family="quasibinomial",
           weights=weights123)')
temp2<-med(data=data.bin,n=2,nonlinear=TRUE)

\donttest{
 #multiple predictors (categorical and continuous predictors)
 x=weight_behavior[,c(3,5:14)]
 pred=weight_behavior[,c(2,4)]
 y=weight_behavior[,15]
 data.b.b.2.3<-data.org(x,y,mediator=4:11,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 temp1.2<-med(data.b.b.2.3,n=2)
 temp2.2<-med(data.b.b.2.3,n=2,nonlinear=TRUE)

 #multivariate responses
 x=weight_behavior[,c(2:3,5:14)]
 pred=weight_behavior[,4]
 y=weight_behavior[,c(1,15)]
 data.b.b.2.4<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 temp1.3<-med(data.b.b.2.4,n=2)
#or use the self defined function
 temp1.3<-med(data.b.b.2.4,n=2,custom.function =c('glm(responseY~.,
              data=dataset123,family="gaussian",weights=weights123)', 
              'glm(responseY~.,data=dataset123,family="quasibinomial",
              weights=weights123)'))
 temp2.3<-med(data.b.b.2.4,n=2,nonlinear=TRUE)

#continuous y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,1]
 data.cont<-data.org(x,y,pred=pred,mediator=5:12,jointm=list(n=1,j1=7:9), 
                     predref="M",alpha=0.4,alpha2=0.4)
 temp3<-med(data=data.cont,n=2) 
 temp4<-med(data=data.cont,n=2,nonlinear=TRUE) 

##continuous x
#binary y
 x=weight_behavior[,3:14]
 pred=weight_behavior[,2]
 y=weight_behavior[,15]
 data.contx<-data.org(x,y,pred=pred,mediator=4:10,alpha=0.4,alpha2=0.4)
 temp5<-med(data=data.contx,n=2)
#or use the self defined function
 temp5<-med(data=data.contx,n=2,custom.function ='glm(responseY~.,
            data=dataset123,family="quasibinomial",weights=weights123)')
 temp6<-med(data=data.contx,n=2,nonlinear=TRUE,nu=0.05)

#continuous y
x=weight_behavior[,3:14]
y=weight_behavior[,1]
pred=weight_behavior[,2]
data.contx<-data.org(x,y,pred=pred,contmed=c(11:12),binmed=c(6,10),
                     binref=c(1,1),catmed=5,catref=1,
                     alpha=0.4,alpha2=0.4)

temp7<-med(data=data.contx,n=2) 
temp8<-med(data=data.contx,n=2,nonlinear=TRUE,nu=0.05) 

##Surv class outcome (survival analysis)
data(cgd1)       #a dataset in the survival package
x=cgd1[,c(4:5,7:12)]
pred=cgd1[,6]
status<-ifelse(is.na(cgd1$etime1),0,1)
y=Surv(cgd1$futime,status)          
#for continuous predictor
data.surv.contx<-data.org(x,y,pred=pred,mediator=1:ncol(x),      
                          alpha=0.5,alpha2=0.5)
temp9.contx<-med(data=data.surv.contx,n=2,type="lp")

#close to mart results when use type="lp"
temp9.contx
temp10.contx<-med(data=data.surv.contx,n=2,nonlinear=TRUE)  
#results in the linear part unit
temp10.contx

#for binary predictor
x=cgd1[,c(5:12)]
pred=cgd1[,4]
data.surv.binx<-data.org(x,y,pred=pred,mediator=1:ncol(x),   
                    alpha=0.4,alpha2=0.4)
temp9.binx<-med(data=data.surv.binx,n=2,type="lp") 
temp9.binx
temp10.binx<-med(data=data.surv.binx,n=2,nonlinear=TRUE)  
temp10.binx}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
