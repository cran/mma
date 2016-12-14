\name{boot.med}
\alias{boot.med}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical Inference on Mediation Analysis with Continuous or Binary Predictor
}
\description{
To make inferences on the mediation effects when the predictor is continuous or binary.
}
\usage{
boot.med(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,
         contm=data$contm,catm=data$catm,jointm=data$jointm,margin=1,
         n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,D=3,
         distn=NULL,family1=NULL,n2=50,weight=rep(1,nrow(x)),
         refy=NULL,x.new=x,binpred=data$binpred,type=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
the list of result from data.org that organize the covariates, mediators, predictor and outcome. If data is FALSE, then need to set x, y, dirx, contm, catm, and jointm.
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
  \item{weight}{
the weight for observations.
}
  \item{refy}{
if y is binary, the reference group of y.
}
  \item{x.new}{
of the same format as x, with a new set of predictor and corresponding covariates on which to calculate the mediation effects.
}
  \item{binpred}{
if TRUE, the predictor is binary.
}
  \item{type}{
the type of prediction when y is class Surv. Is "risk" if not specified.
}
}
\details{
The mediators are not tested in this function.  data.org should be used for the tests and data organizing, and then the resulting list from data.org can be used directly to define the arguments in this function. boot.med considers all variables in x as mediators or covariates in the final model and all variables identified by contm, binm, catm, or jointm as mediators.
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
Qingzhao Yu \email{qyu@lsuhsc.edu}}

\seealso{
\code{"\link[=med]{med}"} just estimate the mediation effects.
}
\examples{
data("weight_behavior")
##binary x
#binary y
x=weight_behavior[,2:14]
y=weight_behavior[,15]
data.bin<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
         binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-boot.med(data=data.bin,n=2,n2=4,seed=1)
temp2<-boot.med(data=data.bin,n=2,n2=4,nu=0.05,seed=1,nonlinear=TRUE)

#continuous y
x=weight_behavior[,2:14]
y=weight_behavior[,1]
data.cont<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
                    binref=c(1,1),catmed=6,catref=1, predref="M",
                    family1=gaussian(link="identity"),alpha=0.4,alpha2=0.4)
\donttest{
temp3<-boot.med(data=data.cont,n=2,n2=4,seed=1) 
temp4<-boot.med(data=data.cont,n=2,n2=4,nu=0.05,seed=1, nonlinear=TRUE)} 

##continuous x
#binary y
x=weight_behavior[,2:14]
y=weight_behavior[,15]
data.contx<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
                     binref=c(1,1),catmed=6,catref=1,alpha=0.4,alpha2=0.4)
\donttest{
temp5<-boot.med(data=data.contx,seed=1,n=1,n2=2)
plot(temp5,vari="exercises",xlim=c(0,30))
temp6<-boot.med(data=data.contx,seed=1,n=1,refy=0,nonlinear=T,n2=2)}

#continuous y
x=weight_behavior[,2:14]
y=weight_behavior[,1]
data.contx<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
                     binref=c(1,1),catmed=6,catref=1,alpha=0.4,alpha2=0.4)
\donttest{
temp7<-boot.med(data=data.contx,seed=1,n=1,n2=2) 
temp8<-boot.med(data=data.contx,nonlinear=T,seed=1,n=1,n2=2) }

##Surv class outcome (survival analysis)
data(cgd0)       #a dataset in the survival package
cgd1<-cgd0
x=cgd1[,c(4:12)]
status<-ifelse(is.na(cgd1$etime1),0,1)
y=Surv(cgd1$futime,status)          
#for continuous predictor
data.surv.contx<-data.org(x,y,pred=3,mediator=(1:ncol(x))[-3],      
                          alpha=0.5,alpha2=0.5)
temp9.contx<-boot.med(data=data.surv.contx,seed=1,n=1,n2=2, type="lp") 
summary(temp9.contx)
\donttest{
temp10.contx<-boot.med(data=data.surv.contx,nonlinear=TRUE,seed=1,n=1,n2=2) 
summary(temp10.contx)}

#for binary predictor
data.surv.binx<-data.org(x,y,pred=1,mediator=(1:ncol(x))[-1],   
                    alpha=0.4,alpha2=0.4)
temp9.binx<-boot.med(data=data.surv.binx,seed=1,n=1,n2=2, type="lp") 
summary(temp9.binx)
\donttest{
temp10.binx<-boot.med(data=data.surv.binx,nonlinear=TRUE,seed=1,n=1,n2=2) 
summary(temp10.binx)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Continuous Predictor }% __ONLY ONE__ keyword per line