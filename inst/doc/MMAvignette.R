## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error=TRUE,
  warning=FALSE
)

## ------------------------------------------------------------------------
library(mma)

## ------------------------------------------------------------------------
data("weight_behavior")
#binary predictor x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,"sex"]
 y=weight_behavior[,"overweigh"]
 data.b.b.1<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="M", alpha=0.4,alpha2=0.4)
 summary(data.b.b.1)

## ---- echo=T, results='hide'---------------------------------------------
 data.b.b.2<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4) 
 summary(data.b.b.2)

## ---- echo=T, results='hide'---------------------------------------------
 cgd1<-cgd0
 status<-ifelse(is.na(cgd1$etime1),0,1)
 y=Surv(cgd1$futime,status) 
 x=cgd1[,c(5:7,9:12)]
 pred=cgd1[,4]
 data.b.surv<-data.org(x,y,pred=pred,mediator=1:ncol(x),   
                    alpha=0.4,alpha2=0.4)
 summary(data.b.surv)

## ------------------------------------------------------------------------
 #multivariate predictor
 x=weight_behavior[,c(2:3,5:14)]
 pred=weight_behavior[,4]
 y=weight_behavior[,15]
 data.mb.b <-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                      pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 summary(data.mb.b)

## ------------------------------------------------------------------------
 #multivariate responses
 x=weight_behavior[,c(2:3,5:14)]
 pred=weight_behavior[,4]
 y=weight_behavior[,c(1,15)]
 data.mb.mb<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                      pred=pred,predref="OTHER", alpha=0.4,alpha2=0.4)
 summary(data.mb.mb)

## ---- echo=T, results='hide',eval=F--------------------------------------
#   med.b.b.2<-med(data=data.b.b.2,n=2,nonlinear=TRUE)

## ---- include=F----------------------------------------------------------
 med.b.b.2<-med(data=data.b.b.2,n=2,nonlinear=TRUE)

## ------------------------------------------------------------------------
 med.b.b.1<-med(data=data.b.b.2,n=2)
 med.b.b.1

## ---- fig.show='hold', fig.height=6, fig.width=4-------------------------
plot(med.b.b.1, data.b.b.2,vari="exercises",xlim=c(0,50))
plot(med.b.b.2, data.b.b.2,vari="sports")

## ---- echo=T, eval=F,results='hide'--------------------------------------
#   med.b.surv.1<-med(data=data.b.surv,n=2,type="lp")
#    #close to mart results when use type="lp"
#   med.b.surv.2<-med(data=data.b.surv,n=2,nonlinear=TRUE)
#   #results in the linear part unit

## ---- include=F----------------------------------------------------------
 med.b.surv.1<-med(data=data.b.surv,n=2,type="lp") 
  #close to mart results when use type="lp"
 med.b.surv.2<-med(data=data.b.surv,n=2,nonlinear=TRUE)  
 #results in the linear part unit

## ---- fig.show='hold', fig.height=5, fig.width=7-------------------------
 bootmed.b.b.1<-boot.med(data=data.b.b.2,n=2,n2=4)
 summary(bootmed.b.b.1)

## ---- fig.show='hold', fig.height=5, fig.width=7, echo=T, eval=F, results='hide'----
#   bootmed.b.b.2<-boot.med(data=data.b.b.2,n=2,n2=4,nu=0.05,nonlinear=TRUE)

## ---- fig.show='hold', fig.height=5, fig.width=7, include=F--------------
 bootmed.b.b.2<-boot.med(data=data.b.b.2,n=2,n2=4,nu=0.05,nonlinear=TRUE)

## ---- fig.show='hold', fig.height=5, fig.width=7-------------------------
 summary(bootmed.b.b.2)

## ---- fig.show='hold', fig.height=7, fig.width=5-------------------------
 plot(bootmed.b.b.1,vari="exercises",xlim=c(0,50))

## ---- fig.show='hold', fig.height=7, fig.width=5-------------------------
 plot(bootmed.b.b.1,vari="sports")

## ---- fig.show='hold', fig.height=5, fig.width=7-------------------------
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 mma.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 summary(mma.b.b.glm)

## ---- fig.show='hold', fig.height=5, fig.width=7, echo=T, eval=F,results='hide'----
#   mma.b.b.mart<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
#                      catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,nonlinear=TRUE,n=2,n2=5)
#   summary(mma.b.b.mart)

## ---- include=F----------------------------------------------------------
 mma.b.b.mart<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,nonlinear=TRUE,n=2,n2=5)

## ---- fig.show='hold', fig.height=5, fig.width=7-------------------------
 summary(mma.b.b.mart)

## ---- fig.height=7, fig.width=5------------------------------------------
plot(mma.b.b.mart,vari="exercises")
plot(mma.b.b.glm,vari="sweat")


