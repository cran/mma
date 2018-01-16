\name{summary.mma}
\alias{summary.mma}
\alias{print.summary.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of an mma project
}
\description{
Compute the estimations, standard deviations and confidence intervals of the mediation effects.
}
\usage{
\method{summary}{mma}(object,..., alpha=0.05, plot=TRUE, RE=FALSE,quant=T) 
\method{print}{summary.mma}(x,...,digit=3) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
a mma object created initially call to mma, boot.met.binx, or boot.met.contx.
}
 \item{x}{
a summary.mma object created initially call to summary.mma.
}
  \item{\dots}{
other arguments passed to the print function.
}
  \item{alpha}{
the alpha level for confidence interval.}
  \item{plot}{
default is TRUE, if ture, draw a barplot of the relative mediation effects with confidence intervals.}
  \item{RE}{
default is FALSE, if ture, show the inferences on relative mediation effects.}
  \item{quant}{
default is TRUE, if ture, draw the confidence intervals of relative effects using quantile.}
  \item{digit}{
the number of decimal digits to keep.
}
}
\details{
summary.mma gives a list of the estimations and summary statistics based on the bootstrap results. If plot=T, draw a barplot of the relative effects of the direct effect of the predictor and indirect effects of the mediators. Relative effects is defined as the (in)direct effect divided by the total effect. The plot is arranged in order from the largest to the smallest relative effect.
}
\value{
The function summary.mma return a list. The first item, results, is the list for mediation effects, and the second item, re, is the list for relative effects. Under them, results have the items ie, te and de; re has the items ie and de. In each of the items, est is the estimation of the corresponding (relative) mediation effects based on the whole data, mean is the average estimated (relative) effects from the bootstrap samples, sd is the standard deviation of the estimates from the bootstrap sample. upbd and lwbd are the upper and lower bound of the confidence interval of the estimation using parametric method from the bootstrap sample, upbd_q and lwbd_q are the corresponding quantiles of the estimation from the bootstrap sample.
  \item{ie }{a matrix of statistics inference on the (relative) indirect effects from the mma object.}
  \item{te }{statistics inference on the total effects from the mma object.}
  \item{de }{statistics inference on the (relative) direct effects from the mma object.}
If plot=T, draw a barplot of the relative mediation effects.
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

\seealso{
\code{"\link[=mma]{mma}"},\code{"\link[=boot.med]{boot.med}"} 
}
\examples{
data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
  temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 summary(temp.b.b.glm)
}

%\keyword{ Summary }
