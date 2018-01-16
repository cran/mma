\name{plot.med}
\alias{plot.med}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the mediation effect on the fitted med object
}
\description{
Plot the marginal effect of the selected variable on the outcome, and the marginal effect of the predictor on the selected variable.
}
\usage{
\method{plot}{med}(x,data,...,vari,xlim=range(data$x[,vari],na.rm=T))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a med_iden object created initially call to med.
}
  \item{x}{
a med object created initially call to med, med.binx, or med.contx.
}
  \item{vari}{
  an indices or the name of the variable to plot.
  }
  \item{xlim}{
  the range of the variable to be plotted.
  }
  \item{\dots}{
other arguments passed to the plot function.
}
}
\details{
plot.med plots the marginal effect of the selected variable on the outcome, and the marginal effect of the predictor on the selected variable. If the predictor is binary, draw a histogram or boxplot of the marginal density of the variable at each different value of the predictor.
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=med]{med}"}
}
\examples{
#library(plotrix)
data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
  binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2,seed=1)
temp2<-med(data=data.bin,n=2,seed=1,nonlinear=TRUE)
plot(temp1,data.bin,vari="exercises",xlim=c(0,50))
plot(temp2,data.bin,vari="sports")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
