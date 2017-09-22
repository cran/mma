\name{data.org}
\alias{data.org}
\title{
Data Organization and Identify Potential Mediators
}
\description{
Do a preliminary data analysis to identify potential mediators and covariates. Each variable listed in jointm is forced in the final estimation model as a mediator. Also organize the data into a format that can be directly used for the mediation analysis functions.
}
\usage{
data.org(x,y,pred,mediator=NULL,contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
                   catref=NULL,jointm=NULL,refy=NULL, family1=NULL,
                   predref=NULL,alpha=0.1,alpha2=0.1,testtype=1, weight=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates
}
  \item{y}{
the vector of outcome variable.The outcome can be binary, continuous, or of "Surv" class (see survival package for help).
}
  \item{pred}{
the column number that locates the predictor in x or the predictor name.
}
  \item{mediator}{
the list of mediators (column numbers in x or by variable names). The mediators to be checked can be identified by "contmed", "binmed" and "catmed", or by this argument, "mediator", where binary and categorical mediators in x are identified as factors or characters, the reference group is the first level of the factor or factorized character. if a mediator has only two unique values, the mediator is identified as binary. If the reference groups need to be changed, the binary or categorical mediators can be listed in binmed or catmed, and the corresponding reference group in binref or catref.
}
  \item{contmed}{
a vector of variable names or column numbers that locate the potential continuous mediators in x.
}
  \item{binmed}{
a vector of column numbers that locate the potential binary mediators in x.
}
  \item{binref}{
the defined reference groups of the binary potential mediators in binmed. The first levels of the mediators if is null.
}
  \item{catmed}{
a vector of variable names or column numbers that locate the potential categorical mediators in x. The first levels of the mediators if is null.
}
  \item{catref}{
the defined reference groups of the categorical potential mediators in catmed.
}
  \item{jointm}{
a list that identifies the mediators that need to be considered jointly, where the first item indicates the number of groups of mediators to be considered jointly, and each of the following items identifies the variable names or column numbers of the mediators in x for each group of joint mediators.
}
  \item{refy}{
if y is binary, the reference group of y. The default is the first level of as.factor(y). y for the reference group is assigned as 0.
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component in generalized linear model.  The default value of family1 is binomial(link = "logit") for binary y, and gaussian(link="identity") for continuous y.
}
  \item{predref}{
if the predictor is binary, identify the reference group of the binary predictor. The default is the first level of the predictor. The value of the predictor is 0 for the reference grop.
}
  \item{alpha}{
the significance level at which to test if the potential mediators (identified by contmed, binmed, and catmed) can be used as a covariate or mediator in estimating y when all variables in x are included in the model.  The default value is alpha=0.1
}
  \item{alpha2}{
the significant level at which to test if a potential mediator is related with the predictor. The default value is alpha2=0.1.
}  
  \item{testtype}{if the testtype is 1 (by default), covariates/mediators are identified using full model; if the testtype is 2, covariates/mediators are tested one by one in models with the predictor only.
}
  \item{weight}{the weight for data analysis, by default is rep(1,length(y)).
}
}

\value{
data.org returns a list with the organized data and identifiers of the potential mediators in the organized data set.
  \item{x }{the organized data frame that include all potential mediators and covariates that should be used to estimate the outcome.}
  \item{dirx }{the column number of x that locates the predictor.}
  \item{contm }{the column numbers of x that locate the potential continuous mediators.}
  \item{binm }{when the predictor is continuous, binm gives the column numbers of x that locate the potential binary mediators.}
  \item{catm }{when the predictor is binary, catm gives the column numbers of x that locate the potential binary or categorical mediators; when the predictor is continuous, catm gives a list where the first item is the number of potential categorical mediators, and the following items give the column numbers of each binarized categorical mediator in x.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators, and each of the following items identifies the column numbers of the mediators in the newly organized x for each group of joint mediators.}
  \item{y }{the vector of outcome.}
  \item{fullmodel}{the full linear model fitted with all potential mediators and covariates.}
  \item{rela}{p-values of tests on the realtionship between the predictor and each potiential mediator.}
  \item{surv}{True if y is of class Surv.}
  \item{P1}{If testtype=2, P1 gives the p-value of the corresponding variables in predicting the outcome when only the variable and predictor are covariates in the model.}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi:10.1037/0022-3514.51.6.1173>. The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations. J. Pers Soc Psychol, 51(6), 1173-1182.
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\note{
All other variables in x but not identified by mediator, contmed, binmed, or catmed are forced in the final model as covariates.  Compared with data.org, joint mediators are considered in this function. Every variable in the jointm should be listed in contmed, binmed, or catmed, and these variables are forced to be included as mediators for further mediation analysis. A variables can be included in more than one groups of joint mediators in jointm.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=data.org]{data.org}"}  that does not consider joint mediators, which can be added freely in the mediation analysis functions later.
}
\examples{
data("weight_behavior")
#binary x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.b.b.2.1<-data.org(x,y,mediator=6:13,jointm=list(n=1,j1=c(6,8,10)),
                        pred=2,predref="M", alpha=0.4,alpha2=0.4)
 summary(data.b.b.2.1)
 #Or you can specify the potential mediators and change the reference 
 #group for binary or categorical mediators. In the following code,
 #potential continuous mediators are columns 8,9,10,12, and 13 of x,
 #binary mediators are columns 7 and 11, and categorical mediator is
 #column 6 of x with 1 to be the reference group for all categorical
 #and binary mediators. 
  data.b.b.2<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,jointm=list(n=1,j1=c(6,8,10)),
   predref="M",alpha=0.4,alpha2=0.4) 
  summary(data.b.b.2)
 #use the mediator argument instead of contmet, binmed and catmed
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.b.c.2<-data.org(x,y,pred=2,mediator=6:13,jointm=list(n=1,j1=8:10), 
   predref="M",alpha=0.4,alpha2=0.4)
 summary(data.b.c.2)
 
#continuous x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.c.b.2<-data.org(x,y,pred=1,mediator=6:13,catref=1,
          jointm=list(n=2,j1=8:10,j2=c(6,8)),alpha=0.4,alpha2=0.4)
 summary(data.c.b.2)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.c.c.2<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,jointm=list(n=2,j1=8:10,j2=c(6,8)),
   alpha=0.4,alpha2=0.4)
 summary(data.c.c.2)
 
#Surv class outcome (survival analysis)
data(cgd0)       #a dataset in the survival package
cgd1<-cgd0
x=cgd1[,c(4:12)]
status<-ifelse(is.na(cgd1$etime1),0,1)
y=Surv(cgd1$futime,status)          
 #for continuous predictor
 #all other variables are considered as potential mediator
 data.surv.contx<-data.org(x,y,pred=3,mediator=(1:ncol(x))[-3],      
                    alpha=0.5,alpha2=0.5) 
 summary(data.surv.contx)

 #for binary predictor
 data.surv.binx<-data.org(x,y,pred=1,mediator=(1:ncol(x))[-1],   
                    alpha=0.4,alpha2=0.4)
 summary(data.surv.binx)
}
\keyword{ Mediator Tests }
