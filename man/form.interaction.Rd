\name{form.interaction}
\alias{form.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create interaction terms of predictor(s) and potential moderator(s).
}
\description{
Create interaction terms of predictor(S) and potential moderator(s). Mainly for linear models.
}
\usage{
form.interaction(x,pred,inter.cov,predref=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
As the x in data.org.
}
  \item{pred}{
The predictor(s) that might have interaction effect with the potential moderators.
  }
  \item{inter.cov}{
The vector of names of potential moderators, which are included in x as covariates.
}
  \item{predref}{
The reference group of the predictor if the predictor is categorical.
  }
}
\details{
form.interaction is used to create interaction terms between predictor(s) and potential moderator(s). The created interaction terms should be used as covariates in linear mediation analysis. The function is not needed for nonlinear mediation method.
}
\author{
Qingzhao Yu \email{qyu@lsuhsc.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=test.moderation]{test.moderation}"}
}
\examples{
data("weight_behavior")
pred=weight_behavior[,3]
x=weight_behavior[,c(2,4:14)]
x=form.interaction(x,pred,inter.cov=c("sports","sweat"),predref=NULL) 
head(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }
