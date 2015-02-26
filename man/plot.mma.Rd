\name{plot.mma}
\alias{plot.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative effects plot of the fitted mma object
}
\description{
Plot the relative effects of the predictor and the mediators in the mma object.
}
\usage{
\method{plot}{mma}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a mma object created initially call to mma, boot.met.binx, or boot.met.contx.
}
  \item{\dots}{
other arguments passed to the plot function.
}
}
\details{
plot.mma plots the relative effects of the direct effect of the predictor and indirect effects of the mediators. Relative effects is defined as the (in)direct effect divided by the total effect. The plot is arranged in order from the largest to the smallest relative effect.
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
 plot(temp.b.b.glm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
