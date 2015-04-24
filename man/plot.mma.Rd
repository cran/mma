\name{plot.mma}
\alias{plot.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative effects plot of the fitted mma object
}
\description{
Plot the marginal effect of the selected variable on the outcome, and the marginal effect of the predictor on the selected variable.
}
\usage{
\method{plot}{mma}(x,...,vari,xlim=range(data$x[,vari],na.rm=T))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a mma object created initially call to mma, boot.met.binx, or boot.met.contx.
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
plot.mma plots the marginal effect of the selected variable on the outcome, and the marginal effect of the predictor on the selected variable. If the predictor is binary, draw a histogram or boxplot of the marginal density of the variable at each different value of the predictor.
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014). "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=mma]{mma}"},\code{"\link[=boot.med.binx]{boot.med.binx}"} , \code{"\link[=boot.med.contx]{boot.med.contx}"}
}
\examples{
data("weight_behavior")
x=weight_behavior[,2:14]
y=weight_behavior[,15]
temp.b.b.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),binref=c(1,1),
   catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, jointm=NULL,margin=1, n=2,
   seed=1,mart=FALSE,nu=0.001,D=3,distn="bernoulli",
   family1=binomial(link = "logit"),n2=2)
plot(temp.b.b.glm,vari="exercises",xlim=c(0,50))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
