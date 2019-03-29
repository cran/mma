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
mma(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
    catref=NULL,jointm=NULL,cova=NULL, refy=rep(NA,ncol(data.frame(y))),
    predref=NULL,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
    nonlinear=F,df1=1,nu=0.001,D=3,distn=NULL,
    family1=as.list(rep(NA,ncol(data.frame(y)))),n2=50,w=rep(1,nrow(x)),
    testtype=1, x.new=NULL, pred.new=NULL, cova.new=NULL, type=NULL,w.new=NULL, 
    all.model=FALSE,xmod=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates. All varaibles in x that are not identified as potential mediators (by mediator, contmed, binmed, catmed, or jointm) are forced in mediation analysis as covariates.
}
  \item{y}{
the vector of outcome variable.
}
  \item{pred}{
the vector/matrix of the predictor(s).
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
  \item{cova}{The data frame of covariates to be used to predict the mediator(s). The covariates in cova cannot be potential mediators. If the covariates are for some specific potential mediators, cova$for.m is the vector of names of potential mediators that use the covariates. Works for continuous predictors (pred) only, also the specified covariates should have no missing if only some mediators uses the covariates. Works mainly for continuous predictors (pred), also the specified covariates should have no missing if only some mediators uses the covariates.
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
  \item{df1}{
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
  \item{w}{
the weight for each observation.
}
  \item{testtype}{if the testtype is 1 (by default), covariates/mediators are identified using full model; if the testtype is 2, covariates/mediators are tested one by one in models with the predictor only.
}
  \item{x.new}{
A new set of predictor and corresponding covariates, of the same format as x (after data.org),  on which to calculate the mediation effects. For continuous predictor only. If is NULL, the mediation effects will be calculated based on the original data set.
}
  \item{pred.new}{
A new set of predictor(s), of the same format as x (after data.org), on which to calculate the mediation effects. For continuous predictor only.
}
  \item{cova.new}{
a new set of covaraite(s).
}
  \item{type}{
the type of prediction when y is class Surv. Is "risk" if not specified.
}
  \item{w.new}{
the weights for new.x.
}
  \item{all.model}{
save all the fitted model from bootstrap samples if TRUE. This needs to be true to make inferences on moderation effects. This argument cannot be used with the function \code{"\link[=mma.par]{mma.par}"}.
}
  \item{xmod}{
If there is a moderator, xmod gives the moderator's name in cova and/or x.
}
}
\details{
mma first tests if the potential mediators defined by binm, contm, and catm should be treated as mediators or covariates (if none, the variable will be deleted from further analysis). All variables identified by jointm are treated as mediators. All other variables in x that are not tested are treated as covariates. Then mma does the mediation effects estimation and inference on the selected variables.  
}
\value{
Returns an mma object.
  \item{estimation }{list the estimation of ie (indirect effect), te  (total effect), and de (direct effect from the predictor) separately.}  \item{bootsresults }{a list where the first item, ie, is a matrix of n2 rows where each column gives the estimated indirect effect from the corresponding mediator (identified by the column name) from the n2 bootstrap samples; the second item, te, is a vector of estimated total effects from the bootstrap sample; and the 3rd item, de, is a vector of estimated direct effect of the predictor from the bootstrap sample.} 
  \item{model }{a list where the first item, MART, is T if mart is fitted for the final model; the second item, Survival, is T if a survival model is fitted; the third item, type, is the type of prediction when a survival model is fitted; the fourth item, model, is the fitted final full model where y is the outcome and all predictor, covariates, and mediators are the explanatory variables; and the fourth item, best.iter is the number of best iterations if MART is used to fit the final model.}
  \item{data }{a list that contains all the used data: x=x, y=y, dirx=dirx, binm=binm, contm=contm, catm=catm, jointm=jointm, binpred=F.}
  \item{boot.detail }{for continuous predictor only: a list that contains the mediation effects on each row of new.pred: new.pred=new.pred, te1, de1, ie1.}
  \item{all_model}{a list with all fitted models from bootstrap samples if all.model is TRUE.This value cannot be created with the function \code{"\link[=mma.par]{mma.par}"}.}
  \item{all_iter}{if all.model is TRUE, a matrix with all fitted best iterations if MART is used. Each row is from one bootstrap sample.  This value cannot be created with the function \code{"\link[=boot.med.par]{boot.med.par}"}.}
    \item{all_boot}{if all.model is TRUE, a matrix with all bootstrap samples. Each row is one bootstrap sample.  This value cannot be created with the function \code{"\link[=boot.med.par]{boot.med.par}"}.}
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
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\seealso{
\code{"\link[=data.org]{data.org}"} is for mediator tests, and \code{"\link[=med]{med}"} , and \code{"\link[=boot.med]{boot.med}"} for mediation analysis and inferences.
}
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

\donttest{
##Surv class outcome (survival analysis)
data(cgd0)       #a dataset in the survival package
x=cgd1[,c(4:5,7:12)]
pred=cgd1[,6]
status<-ifelse(is.na(cgd1$etime1),0,1)
y=Surv(cgd1$futime,status)          
#for continuous predictor
temp.cox.contx<-mma(x,y,pred=pred,mediator=1:ncol(x),      
                    alpha=0.5,alpha2=0.5,seed=1,n=1,n2=2,type="lp") 
summary(temp.cox.contx)
temp.surv.mart.contx<-mma(x,y,pred=pred,mediator=1:ncol(x),      
                    alpha=0.5,alpha2=0.5,nonlinear=TRUE,seed=1,n2=2) 
summary(temp.surv.mart.contx,ball.use=F)
plot(temp.surv.mart.contx,vari="steroids")
plot(temp.cox.contx,vari="steroids")

#for binary predictor
x=cgd1[,c(5:12)]
pred=cgd1[,4]
temp.cox.binx<-mma(x,y,pred=pred,mediator=1:ncol(x),   
                    alpha=0.4,alpha2=0.4,seed=1,n=1,n2=2,type="lp") 
summary(temp.cox.binx)
temp.surv.mart.binx<-mma(x,y,pred=pred,mediator=1:ncol(x),   
                    alpha=0.4,alpha2=0.4,nonlinear=TRUE,n=1,n2=2) 
summary(temp.surv.mart.binx, RE=T,ball.use=F)
plot(temp.surv.mart.binx,vari="hos.cat")
plot(temp.cox.binx,vari="hos.cat")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Mediator Tests }% __ONLY ONE__ keyword per line
