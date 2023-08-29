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
\method{summary}{mma}(object,..., alpha=0.05, plot=TRUE, RE=FALSE,
quant=FALSE,ball.use=FALSE,bymed=FALSE,win1=0.2) 
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
default is TRUE, if ture, draw a barplot of the mediation effects with confidence intervals.}
  \item{RE}{
default is FALSE, if ture, show the inferences on relative mediation effects.}
  \item{quant}{
default is TRUE, if ture and ball.use is FALSE, draw the confidence intervals of relative effects using quantile.}
  \item{ball.use}{
default is TRUE, if ture, draw the confidence intervals of relative effects using the confidence ball. If both quant and ball.use are false, draw the confidence intervals based on the standard deviaitons from bootstrap estimates.}
  \item{digit}{
the number of decimal digits to keep.
}
  \item{bymed}{
if true, show results by each mediator.
}
  \item{win1}{
The rate of winsarize.
}
}
\details{
summary.mma gives a list of the estimations and summary statistics based on the bootstrap results. If plot=T, draw a barplot of the relative effects of the direct effect of the predictor and indirect effects of the mediators. Relative effects is defined as the (in)direct effect divided by the total effect. The plot is arranged in order from the largest to the smallest relative effect.
}
\value{
The function summary.mma return a list. The first item, results, is the list for mediation effects, and the second item, re, is the list for relative effects. Under them, results have the items ie, te and de; re has the items ie and de. In each of the items, est is the estimation of the corresponding (relative) mediation effects based on the whole data, mean is the average estimated (relative) effects from the bootstrap samples, sd is the standard deviation of the estimates from the bootstrap sample. upbd and lwbd are the upper and lower bound of the confidence interval of the estimation using parametric method from the bootstrap sample, upbd_q and lwbd_q are the corresponding quantiles of the estimation from the bootstrap sample.
  \item{ie }{a matrix of statistics inference on the (relative) indirect effects from the mma object.  est is the estimate using the full sample. mean is the estimate that average over the bootstrap estimates. sd is the standard deviation of the bootstrap estimates. upbd and lwbd are the confidence bounds based on sd. upbd_q and lwbd_q are the confidence bounds based on quantiles of the bootstrap estimates. upbd_b and lwbd_b are condidence ball bounds based on bootstap estimates.}
  \item{te }{statistics inference on the total effects from the mma object.}
  \item{de }{statistics inference on the (relative) direct effects from the mma object.}
If plot=T, draw a barplot of the relative mediation effects.
}
\references{
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
\code{"\link[=mma]{mma}"},\code{"\link[=boot.med]{boot.med}"} 
}
\examples{
data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
  temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 summary(temp.b.b.glm, RE=TRUE, ball.use=FALSE)
 summary(temp.b.b.glm, ball.use=FALSE)
}

%\keyword{ Summary }
