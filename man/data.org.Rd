\name{data.org}
\alias{data.org}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Organization
}
\description{
Do a preliminary data analysis to identify potential mediators and covariates.  Also organize the data into a format that can be directly used for the mediation analysis functions.
}
\usage{
data.org(x, y, pred, contmed = NULL, binmed = NULL, binref = NULL, 
  catmed = NULL, catref = NULL, biny = T, family1 = binomial(link = "logit"), 
  binpred = T, predref = 1, alpha = 0.1, alpha2 = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates
}
  \item{y}{
the vector of outcome variable.
}
  \item{pred}{
the column number that locates the predictor in x.
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
  \item{biny}{
if TRUE, the outcome variable is binary.  The default value of biny is TRUE.
}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component in generalized linear model.  The default value of family1 is binomial(link = "logit"). If y is continuous, family1 can be defined as gaussian(link="identity").
}
  \item{binpred}{
if TRUE, the predictor is binary.  The default value of binpred is TRUE.
}
  \item{predref}{
if binpred is TRUE, identify the reference group of the binary predictor.
}
  \item{alpha}{
the significance level at which to test if the potential mediators (identified by contmed, binmed, and catmed) can be used as a covariate or mediator in estimating y when all variables in x are included in the model.  The default value is alpha=0.1
}
  \item{alpha2}{
the significant level at which to test if a potential mediator is related with the predictor. The default value is alpha2=0.1.
}
}
\value{
 data.org returns a list with the organized data and identifiers of the potential mediators in the organized data set.
  \item{x }{the organized data frame that include the predictor, all potential mediators and covariates that should be used to estimate the outcome.}
  \item{dirx }{the column number of x that locates the predictor.}
  \item{contm }{the column numbers of x that locate the potential continuous mediators.}
  \item{binm }{when the predictor is continuous, binm gives the column numbers of x that locate the potential binary mediators.}
  \item{catm }{when the predictor is binary, catm gives the column numbers of x that locate the potential binary or categorical mediators; when the predictor is continuous, catm gives a list where the first item is the number of potential categorical mediators, and the following items give the column numbers of each binarized categorical mediator in x.}
}
\references{
Baron, R.M., and Kenny, D.A. (1986). The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations. J. Pers Soc Psychol, 51(6), 1173-1182.
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\note{
All other variables in x but not identified by contmed, binmed, or catmed are forced in the final model as covariates.  Compared with data.org2, joint mediators are not considered in this function but can be added freely in the mediation functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=data.org2]{data.org2}"} to consider joint mediators, of which each variable is forced in the final estimation model as a mediator.
}
\examples{
data("weight_behavior")
#binary x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.b.b<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.b.c<-data.org(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,biny=FALSE, 
   family1=gaussian(link="identity"),predref="M",alpha=0.4,alpha2=0.4)
 
#continuous x
 #binary y
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 data.c.b<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,binpred=FALSE,alpha=0.4,alpha2=0.4)
 
 #continuous y
 x=weight_behavior[,2:14]
 y=weight_behavior[,1]
 data.c.c<-data.org(x,y,pred=1,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,binpred=FALSE,
   family1=gaussian(link="identity"),alpha=0.4,alpha2=0.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediator Tests }
