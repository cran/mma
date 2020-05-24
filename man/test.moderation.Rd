\name{test.moderation}
\alias{test.moderation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for moderation effects.
}
\description{
Used on a med object, test if there are moderation effects.
}
\usage{
test.moderation(med1,vari,j=1,kx=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{med1}{
a med object created initially call to med.
}
  \item{vari}{
If med1 was from generalized linear method, vari is the vector of names of potential moderation effects -- the interaction terms formed by the potential moderator(s) times the predictor(s). Those interaction terms have to be included in mediation analysis as covariates.

If med1 was formed by nonlinear method, vari is the vector of names of potential moderators, which have been included in mediation analysis as covariates.
}
  \item{j}{
To check the moderation effect on the jth response.
  }
  \item{kx}{
The moderation effect is with the kx-th predictor(s). kx can be a vector. If kx is null, the moderation effect is with each predictor.
  }
}
\details{
test.moderation is used to test whether direct moderation effect is significant based on the med object from med function (mediation analysis). If geralized linear models were used for mediation anlysis, the moderation effects (formed by the interactions of potential moderators and predictor) should be included in the mediation analysis as covariates. The function test.modeartion will give the significance levels of these interaction terms. If nonlinear models were used for mediation analysis, the potential moderators (but not the interaction terms) should have been included in mediation analysis as covarites. The function test.moderation will give not only significance levels of their interaction with the predictor(s) in generalized linear model, but also the H-statistics by Friedman and Popescue (2008).
}
\value{No return value, called for test results.}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189.

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23.

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11.

Yu, Q., Wu, X., Li, B., and Scribner, R. (2018). <doi:10.1002/sim.7977>. "Multiple Mediation Analysis with Survival Outcomes – With an Application to Explore Racial Disparity in Breast Cancer Survival," Statistics in Medicine.

Yu, Q., Medeiros, KL, Wu, X., and Jensen, R. (2018). <doi:10.1007/s11336-018-9612-2>. "Explore Ethnic Disparities in Anxiety and Depression Among Cancer Survivors Using Nonlinear Mediation Analysis," Psychometrika, 83(4), 991-1006.

Friedman, J.H. and Popescu B.E. (2008) "PREDICTIVE LEARNING VIA RULE ENSEMBLES," Annals of Applied Statistics, 2(3): 916-954.
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=med]{med}"},, \code{"\link[=moderate]{moderate}"
}}
\examples{
data("weight_behavior")
x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,15]
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp2<-med(data=data.bin,n=2,nonlinear=TRUE)
test.moderation(temp2,c("sports","sweat"),j=1,kx=NULL)

%For linear mma, first create the interaction terms
x=cbind(x,form.interaction(x,pred,inter.cov=c("sports","sweat"),predref=NULL))

data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2)
test.moderation(temp1,c("sports","sweat"),j=1,kx=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
