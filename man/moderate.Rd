\name{moderate}
\alias{moderate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate and plot the direct effect of the selected exposure variable at each level of the moderator.
}
\description{
Calculate and plot the direct effect of the selected exposure variable at each level of the moderator.
}
\usage{
moderate(med1,vari,j=1,kx=1,continuous.resolution=100,plot=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{med1}{
The med object from the \link[=med]{med} function.
}
  \item{vari}{
The name of the moderator.
  }
  \item{j}{
The jth response if the response is multiple.
}
  \item{kx}{
The moderate effect is with the kx-th predictor(s). 
  }
  \item{continuous.resolution}{
The number of equally space points at which to evaluate continuous predictors. 
  }
  \item{plot}{
Plot the direct effect at each level of the moderator if ture. 
  }
}
\value{
The \link[=moderate]{moderate} returns a list where the item result is a data frame with two or three elements
  \item{moderator }{the moderator levels.}
  \item{x }{the level of the exposure variable -- available only for continuous exposure and moderate with nonlinear method.}
  \item{de }{the direct effect at the corresonding moderator (and exposure) level(s).}
}
\details{
Calculate and plot the direct effect of the selected exposure variable at each level of the moderator base on the result from the \link[=med]{med} function.
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=form.interaction]{form.interaction}"},
\code{"\link[=test.moderation]{test.moderation}"}
}
\examples{
\donttest{
#nonlinear model
 data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
     binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp2<-med(data=data.bin,n=2,seed=1,nonlinear=TRUE)
 result1=moderate(temp2,vari="race")
 result2=moderate(temp2,vari="age")

#linear model
 data("weight_behavior")
 pred=weight_behavior[,3]
 x=weight_behavior[,c(2,4:14)]
 inter=form.interaction(x,pred,inter.cov=c("race","age"),predref="M") 
 x=cbind(x,inter)
 head(x)
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
     binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-med(data=data.bin,n=2,seed=1)
 result3=moderate(temp1,vari="race")
 result4=moderate(temp1,vari="age")

##with a transformation of continuous moderator
x=weight_behavior[,c(2,4:14)]
x=cbind(x,age2=x[,"age"]^2)
inter=form.interaction(x,pred,inter.cov=c("age","age2"),predref="M") 
x=cbind(x,inter)
head(x)
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2,seed=1)
result5=moderate(temp1,vari="age")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
