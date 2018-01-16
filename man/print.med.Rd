\name{print.med}
\alias{print.med}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print an med object
}
\description{
Print the estimation of mediation effects from an med object: from functions med.
}
\usage{
\method{print}{med}(x,...,digit=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a med object created initially call to med.
}
  \item{\dots}{
other arguments passed to the print function.
}
  \item{digit}{
the number of digits to keep at printing.
}
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
data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                    binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-med(data=data.bin,n=2,seed=1)
 temp2<-med(data=data.bin,n=2,seed=1,nonlinear=TRUE)
 temp1
 print(temp2,digit=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Print }
