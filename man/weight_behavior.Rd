\name{weight_behavior}
\alias{weight_behavior}
\docType{data}
\title{
Weight_Behavior Data Set
}
\description{
This database was obtained from the Louisiana State University Health Sciences Center, New Orleans, by Dr. Richard Scribner.  He explored the relationship between BMI and kids behavior through a survey at children, teachers and parents in Grenada in 2014.  This data set includes 691 observations and 15 variables. 
}
\usage{data(weight_behavior)}
\format{
The data set contains the following variables:

  bmi - body mass index, calculated by weight(kg)/height(cm)^2, numeric
  
  age - children's age in years at the time of survey, numeric
  
  sex - sex of the children, factor
  
  race - race of the children, factor
  
  numpeople - number of people in family, numeric
  
  car - the number of cars in family, numeric
  
  gotosch -  the method used to go to school, factor
  
  snack - eat snack or not in a day, binary
  
  tvhours - number of hours watching TV per week, numeric
  
  cmpthours - number of hours using computer per week, numeric
  
  cellhours - number of hours playing with cell phones per week, numeric
  
  sports - join in a sport team or not, 1: yes; and 2: no
  
  exercises - number of hours of exercises per week, numeric
  
  sweat - number of hours of sweating activities per week, numeric
  
  overweigh - the child is overweighed or not, binary
}
\examples{
data(weight_behavior)
names(weight_behavior)
}
\keyword{Datasets}
