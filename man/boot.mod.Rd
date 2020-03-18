\name{boot.mod}
\alias{boot.mod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical Inference on Mediation Analysis with Continuous or Binary Predictor at different level of the moderator
}
\description{
To make inferences on the mediation effects when the predictor is continuous or binary at different level of the moderator.
}
\usage{
boot.mod(mma1,vari,continuous.resolution=10,
         w=NULL,n=20,
         x.new=NULL,w.new=NULL,pred.new=NULL,cova.new=NULL,xj=1,margin=1,
         xmod=vari,df1=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{mma1}{
the mma project generated from the function mma or boot.med.
}
   \item{vari}{
The name of the moderator. The moderator should be included in the original data as a covariate for the response(in x) and/or for mediators (in cova). 
}
  \item{continuous.resolution}{
The number of equally space points (for continuous moderator) or levels (for categorical moderator) at which to evaluate mediation effects separately.
}
  \item{w}{
the weight for observations.
}
  \item{n}{
the time of resampling in calculating the indirect effects, default is n=20, see Yu et al. (2014).
}
%  \item{seed}{
%set seed to make the calculation repeatable. The default value of seed is from sample(1:1000,1).
%}
  \item{x.new}{
of the same format as x, with a new set of covariates and mediators on which to calculate the mediation effects.
}
  \item{w.new}{
the weights for new.x.
}
  \item{pred.new}{
a new set of predictor(s).
}
  \item{cova.new}{
a new set of covariate for mediator(s).
}
  \item{xj}{
the inference of mediation effects at different levels of moderator is made on the xjth predictor.
}
  \item{margin}{
the change in predictor when calculating the mediation effects, see Yu et al. (2014).
}
  \item{xmod}{
If there is a moderator, xmod gives the moderator's name in cova and/or x.
}
  \item{df1}{
if nonlinear is TURE, natural cubic spline will be used to fit the relationship between the predictor and each mediator.  The df is the degree of freedom in the ns() function, the default is 1.
}
}
\details{
calculate the mediation effects from the xjth predictor to the response variable(s) at each level of the mediator.
}
\value{
Returns an mma.mod object, which is similar to the mma object. Instead of calculating mediation effects for each predictor, the function returns the mediation effects of the xjth predictor at each level of the moderator.
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Yu, Q., Wu, X., Li, B., and Scribner, R. (2018). <doi:10.1002/sim.7977>. "Multiple Mediation Analysis with Survival Outcomes – With an Application to Explore Racial Disparity in Breast Cancer Survival," Statistics in Medicine.

Yu, Q., Medeiros, KL, Wu, X., and Jensen, R. (2018). <doi:10.1007/s11336-018-9612-2>. "Explore Ethnic Disparities in Anxiety and Depression Among Cancer Survivors Using Nonlinear Mediation Analysis," Psychometrika, 83(4), 991-1006. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}}

\seealso{
\code{"\link[=med]{med}"} just estimate the mediation effects.
}
\examples{
\donttest{
#binary x and categorical moderator
#linear model
data("weight_behavior")
pred=weight_behavior[,3]
x=weight_behavior[,c(2,4:14)]
inter=form.interaction(x,pred,inter.cov=c("race"),predref="M") 
x=cbind(x,inter)
y=weight_behavior[,15]
head(x)
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-boot.med(data=data.bin,n=2,n2=4,all.model=T)
temp1.mode=boot.mod(temp1,vari="race",continuous.resolution=100)
#summary(temp1.mode,ball.use = F)
plot2.mma(temp1.mode,vari="exercises",moderator="race")
plot2.mma(temp1.mode,vari="sports",moderator="race")

#nonlinear model
x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,15]
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp2<-boot.med(data=data.bin,n=2,n2=4,nu=0.05,nonlinear=TRUE,all.model=T)
temp2.mode=boot.mod(temp2,vari="race",continuous.resolution=100)
#summary(temp2.mode,ball.use = F)
plot2.mma(temp2.mode,vari="exercises",moderator="race")
plot2.mma(temp2.mode,vari="sports",moderator="race")

#binary x and continuou moderator
x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,1]
inter=form.interaction(x,pred,inter.cov=c("age"),predref="M") 
x=cbind(x,inter)
head(x)
data.cont<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp3<-boot.med(data=data.cont,n=2,n2=4,all.model=T) 
temp3.mode=boot.mod(temp3,vari="age",continuous.resolution=2)
summary(temp3.mode)

plot2.mma(temp3.mode,vari="exercises",moderator="age")
plot2.mma(temp3.mode,vari="sports",moderator="age")

x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,1]
data.cont<-data.org(x,y,pred=pred,mediator=5:12,jointm=list(n=1,j1=7:9), 
                     predref="M",alpha=0.4,alpha2=0.4)
temp4<-boot.med(data=data.cont,n=2,n2=4,nu=0.05, nonlinear=TRUE,all.model=T)
temp4.mode=boot.mod(temp4,vari="age",continuous.resolution=2,n=20)
summary(temp4.mode)
summary(temp4.mode,bymed=T)

#categorical m and continuous x
x=weight_behavior[,3:14]
pred=weight_behavior[,2]
inter=form.interaction(x,pred,inter.cov=c("race")) 
cova=data.frame(inter)
x=cbind(x,cova)
head(x)
y=weight_behavior[,15]
data.contx<-data.org(x,y,pred=pred,mediator=4:10,alpha=0.4,alpha2=0.4,cova=cova)
temp5<-boot.med(data=data.contx,n=1,n2=2,cova=cova,all.model=T)
temp5.mode=boot.mod(temp5,vari="race",continuous.resolution=10)
summary(temp5.mode)

x=weight_behavior[,3:14]
pred=weight_behavior[,2]
y=weight_behavior[,15]
data.contx<-data.org(x,y,pred=pred,mediator=4:10,alpha=0.4,alpha2=0.4,cova=cova)
temp6<-boot.med(data=data.contx,n=1,refy=0,nonlinear=T,n2=3,
                cova=data.frame(race=x$race),all.model = T)
temp6.mode=boot.mod(temp6,vari="race",continuous.resolution=10)
summary(temp6.mode)

#plot2.mma(temp5.mode,vari="car",moderator="race")
#plot2.mma(temp6.mode,vari="car",moderator="race")

#continuous m and x
x=weight_behavior[,3:14]
pred=weight_behavior[,2]
inter=form.interaction(x,pred,inter.cov=c("numpeople")) 
x=cbind(x,inter)
head(x)
y=weight_behavior[,15]
data.contx<-data.org(x,y,pred=pred,mediator=4:10,alpha=0.4,alpha2=0.4)
temp5.1<-boot.med(data=data.contx,n=1,n2=2,all.model=T)
temp5.mode1=boot.mod(temp5.1,vari="numpeople",continuous.resolution=2)
summary(temp5.mode1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Moderation Effect with Mediation Analysis }
%\keyword{ Continuous Predictor }% __ONLY ONE__ keyword per line
