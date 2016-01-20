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
Yu, Q., Fan, Y., and Wu, X. (2014). "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

\seealso{
\code{"\link[=mma]{mma}"},\code{"\link[=boot.med.binx]{boot.med.binx}"} , \code{"\link[=boot.med.contx]{boot.met.contx}"}
}
\examples{
data("weight_behavior")
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),mart=FALSE,nu=0.001,
   D=3,distn="bernoulli",family1=binomial(link = "logit"),n2=2)
 summary(temp.b.b.glm)
}

%\keyword{ Summary }
