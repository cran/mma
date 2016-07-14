\name{print.mma}
\alias{print.mma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a mma object
}
\description{
Print the estimation of mediation effects from an mma object.
}
\usage{
\method{print}{mma}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a mma object created initially call to mma, boot.med.binx, or boot.med.contx.
}
  \item{\dots}{
other arguments passed to the print function.
}
}
\value{
Return a LIST
  \item{ie }{estimation of the indirect effects from the mma object.}
  \item{te }{estimation of the total effect from the mma object.}
  \item{de }{estimation of the direct effect from the mma object.}
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi: 10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189.
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=mma]{mma}"},\code{"\link[=boot.med]{boot.med}"}
}
\examples{
data("weight_behavior")
 x=weight_behavior[,2:14]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=2,contmed=c(8:10,12:13),binmed=c(7,11),
   binref=c(1,1),catmed=6,catref=1,predref="M",alpha=0.4,alpha2=0.4, 
   jointm=NULL,margin=1, n=2,seed=sample(1:1000,1),nonlinear=FALSE,nu=0.001,
   D=3,n2=2)
 print(temp.b.b.glm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Print }
