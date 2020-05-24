\name{summary.med_iden}
\alias{summary.med_iden}
\alias{print.summary.med_iden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for class "med_iden".
}
\description{
Compute the estimations, standard deviations and confidence intervals of the mediation effects.
}
\usage{
\method{summary}{med_iden}(object,...,only=FALSE) 
\method{print}{summary.med_iden}(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
a med_iden object created initially call to data.org.
}
 \item{x}{
a summary.med_iden object created initially call to summary.med_iden
}
 \item{\dots}{
other arguments passed to the print function.
}
 \item{only}{
if only=T, show test results for selected covariates and mediators only.
}
}
\details{
summary.med_iden gives a list of identified mediators, covariates and the test results. 
}
\value{
The function summary.med_iden return a list of covariates and mediators as identified by a series tests.
  \item{mediator }{variable names of the identified mediators, either siginicant in both full model and in relate to the predictor, or being a member of the pre-identified joint mediators.}
  \item{covariate }{variable names of covariates: being significant in the full model but not significantly relate with the predictor.}
  \item{tests }{statistical test results.}
  \item{results }{the original object.}
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
 data.b.b.2<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="M", alpha=0.4,alpha2=0.4)
 summary(data.b.b.2)
}

%\keyword{ Summary }
