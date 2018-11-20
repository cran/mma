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
test.moderation is used to test whether certain moderation effect is significant based on the med object from med function (mediation analysis). If geralized linear models were used for mediation anlysis, the moderation effects (formed by the interactions of potential moderators and predictor) should be included in the mediation analysis as covariates. The function test.modeartion will give the significance levels of these interaction terms. If nonlinear models were used for mediation analysis, the potential moderators (but not the interaction terms) should have been included in mediation analysis as covarites. The function test.moderation will give not only significance levels of their interaction with the predictor(s) in generalized linear model, but also the H-statistics by Friedman and Popescue (2008).
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189.

Friedman, J.H. and Popescu B.E. (2008) "PREDICTIVE LEARNING VIA RULE ENSEMBLES," Annals of Applied Statistics, 2(3): 916-954.
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=med]{med}"}
}
\examples{
data("weight_behavior")
x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,15]
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp2<-med(data=data.bin,n=2,seed=1,nonlinear=TRUE)
test.moderation(temp2,c("sports","sweat"),j=1,kx=NULL)

%For linear mma, first create the interaction terms
x=form.interaction(x,pred,inter.cov=c("sports","sweat"),predref=NULL) 
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2,seed=1)
test.moderation(temp1,c("sports","sweat"),j=1,kx=NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
