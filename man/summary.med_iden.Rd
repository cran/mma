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
\method{summary}{med_iden}(object,...) 
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
