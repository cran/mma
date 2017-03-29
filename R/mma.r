#to organize data
data.org<-function(x,y,pred,mediator=NULL,contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
                   catref=NULL,jointm=NULL,refy=NULL, family1=NULL,
                   predref=NULL,alpha=0.1,alpha2=0.1,testtype=1)
{cattobin<-function(x,cat1,cat2=rep(1,length(cat1)))
{ ad1<-function(vec)
{vec1<-vec[-1]
 vec1[vec[1]]<-1
 vec1
}
 dim1<-dim(x)
 catm<-list(n=length(cat1))
 g<-dim1[2]-length(cat1)
 ntemp<-names(x)[cat1]
 j<-1
 for (i in cat1)
 {a<-x[,i]
  d<-rep(0,dim1[1])
  b<-sort(unique(a[a!=cat2[j]]))
  l<-1
  for (k in b)
  {d[a==k]<-l
   l<-l+1}
  d[a==cat2[j]]<-l
  f<-matrix(0,dim1[1],l-1) 
  colnames(f)<-paste(ntemp[j],b,sep=".") #changed for error info
  hi<-d[d!=l & !is.na(d)]
  f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
  f[is.na(d),]<-NA
  x<-cbind(x,f)
  catm<-append(catm,list((g+1):(g+l-1)))
  g<-g+length(b)
  j<-j+1
 }
 x<-x[,-cat1]
 list(x=x,catm=catm)
}
colnum<-function(vec,cutx) 
{z<-vec
 for (i in 1:length(vec))
   z[i]<-vec[i]-sum(vec[i]>cutx)
 z}

surv<-F
if(class(y)=="Surv")
  {surv<-T
   biny<-F}
else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
{biny=T
if(is.null(family1))
  family1 = binomial("logit")
#if(is.null(distn))
#  distn = "bernoulli"
if(!is.null(refy))
  y<-ifelse(y==refy,0,1)
else
  y<-ifelse(as.factor(y)==levels(as.factor(y))[1],0,1)
}
else
{biny=F
if(is.null(family1))
  family1 = gaussian(link = "identity")
#if(is.null(distn))
#  distn = "gaussian"
}

xnames<-colnames(x)
if(is.character(pred))
  pred<-grep(pred,xnames)

if(is.character(x[,pred]) | is.factor(x[,pred]) | nlevels(as.factor(x[,pred]))==2)
 {binpred=T
  if(!is.null(predref))
      x[,pred]<-as.factor(ifelse(x[,pred]==predref,0,1))
  else
      {x[,pred]<-as.factor(x[,pred])
       x[,pred]<-as.factor(ifelse(x[,pred]==levels(x[,pred])[1],0,1))}
 }
else
  binpred=F

if(is.character(contmed))
  contmed<-unlist(sapply(contmed,grep,xnames))
if(is.character(binmed))
  binmed<-unlist(sapply(binmed,grep,xnames))
if(is.character(catmed))
  catmed<-unlist(sapply(catmed,grep,xnames))
if(!is.null(jointm))
  for (i in 2:length(jointm))
   if(is.character(jointm[[i]]))
     jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))

if(!is.null(binmed) & is.null(binref))
  for(i in binmed)
  {x[,i]=as.factor(x[,i])
   binref=c(binref,levels(x[,i])[1])}

if(!is.null(catmed) & is.null(catref))
  for(i in catmed)
  {x[,i]=as.factor(x[,i])
   catref=c(catref,levels(x[,i])[1])}

if(!is.null(mediator))   #add a mediator argument
 {if(is.character(mediator))
   mediator<-unlist(sapply(mediator,grep,xnames))
  for (i in 1:length(mediator))
    {if(is.character(x[,mediator[i]]))
      x[,mediator[i]]<-as.factor(x[,mediator[i]])
     if(is.factor(x[,mediator[i]]) | nlevels(as.factor(x[,mediator[i]]))==2)
       if(nlevels(as.factor(x[,mediator[i]]))==2)
       {if(sum(binmed==mediator[i])==0)
        {x[,mediator[i]]<-as.factor(x[,mediator[i]])
         binmed<-c(binmed,mediator[i])
         binref<-c(binref,(levels(x[,mediator[i]]))[1])
         }}
      else
       {if(sum(catmed==mediator[i])==0)
        {catmed<-c(catmed,mediator[i])
         catref<-c(catref,levels(x[,mediator[i]])[1])}}
     else if(sum(contmed==mediator[i])==0 & sum(catmed==mediator[i])==0)
        contmed<-c(contmed,mediator[i])}
 }

if(!is.null(jointm))  #identify mediators that should be jointly considered
{joint<-NULL 
 for (i in 2:(jointm[[1]]+1))
   joint=c(joint,jointm[[i]])
 joint1<-unique(joint)
 
 if(!is.null(contmed))
 {cont1<-rep(F,length(contmed))
  for (i in 1:length(contmed))
    cont1[i]<-ifelse(sum(contmed[i]==joint1)>0,T,F)
 }
 if(!is.null(binmed))
 {bin1<-rep(F,length(binmed))
  for (i in 1:length(binmed))
    bin1[i]<-ifelse(sum(binmed[i]==joint1)>0,T,F)
 }
 if(!is.null(catmed))
 {cat1<-rep(F,length(catmed))
  for (i in 1:length(catmed))
    cat1[i]<-ifelse(sum(catmed[i]==joint1)>0,T,F)
 }
}
else
{if(!is.null(contmed))
  cont1<-rep(F,length(contmed))
 
 if(!is.null(binmed))
   bin1<-rep(F,length(binmed))
 
 if(!is.null(catmed))
   cat1<-rep(F,length(catmed))
}

if(!is.null(binmed))
{j<-1
for (i in binmed)
{x[,i]<-ifelse(x[,i]==binref[j],0,1)
 j<-j+1}
}

if(!is.null(catmed))
{tempx<-cattobin(x,catmed,catref)
 newx<-tempx$x
 catnewx<-tempx$catm}
else newx<-x
#delete variables that are not significant
if(surv)
  fullmodel<-summary(coxph(y~.,data=data.frame(newx)))
else
  fullmodel<-summary(glm(y~.,data=data.frame(newx),family=family1))
xname<-names(x)
if(surv)
  fullmodel1<-coxph(y~.,data=data.frame(x))
else
  fullmodel1<-glm(y~.,data=data.frame(x),family=family1) #use type three error to test the full model
type3<-Anova(fullmodel1,type="III")
xnames3<-rownames(type3)
P1<-matrix(NA,length(xnames3),1)
rownames(P1)<-xnames3
if(surv)
  temp.fullmodel1<-coxph(y~x[,pred])
else
  temp.fullmodel1<-glm(y~x[,pred],family=family1) #use type three error to test the full model

P1[grep(xname[pred],xnames3),]<-Anova(temp.fullmodel1,type="III")[1,3]

covr.cont<-rep(F,length(contmed))
if(!is.null(contmed))
{for (i in 1:length(contmed))
  {if(testtype==1)
    covr.cont[i]<-ifelse(type3[xnames3==xname[contmed[i]],3]<alpha,T,F)
   else if(testtype==2)
   {temp.data<-x[,c(contmed[i],pred)]
    if(surv)
     temp.fullmodel1<-coxph(y~.,data=data.frame(temp.data))
    else
     temp.fullmodel1<-glm(y~.,data=data.frame(temp.data),family=family1) #use type three error to test the full model
    temp.type3<-Anova(temp.fullmodel1,type="III")
    temp.p1<-temp.type3[rownames(temp.type3)==xname[contmed[i]],3]
    P1[grep(xname[contmed[i]],xnames3),]<-temp.p1
    covr.cont[i]<-ifelse(temp.p1<alpha,T,F)
   }
  }
  #covr.cont[i]<-ifelse(fullmodel$coef[dimnames(fullmodel$coef)[[1]]==xname[contmed[i]],4]<alpha,T,F)
 covr.cont<-ifelse(covr.cont|cont1,T,F)
 cont2<-cont1[covr.cont]
 contmed1<-contmed[covr.cont]
}

covr.bin<-rep(F,length(binmed))
if(!is.null(binmed))
{for (i in 1:length(binmed))
 {if(testtype==1)
  covr.bin[i]<-ifelse(type3[xnames3==xname[binmed[i]],3]<alpha,T,F)
 else if(testtype==2)
 {temp.data<-x[,c(binmed[i],pred)]
  if(surv)
   temp.fullmodel1<-coxph(y~.,data=data.frame(temp.data))
  else
   temp.fullmodel1<-glm(y~.,data=data.frame(temp.data),family=family1) #use type three error to test the full model
 temp.type3<-Anova(temp.fullmodel1,type="III")
 temp.p1<-temp.type3[rownames(temp.type3)==xname[binmed[i]],3]
 P1[grep(xname[binmed[i]],xnames3),]<-temp.p1
 covr.bin[i]<-ifelse(temp.p1<alpha,T,F)
  }
 }
  #covr.bin[i]<-ifelse(fullmodel$coef[dimnames(fullmodel$coef)[[1]]==xname[binmed[i]],4]<alpha,T,F)
 covr.bin<-ifelse(covr.bin+bin1>0,T,F) 
 bin2<-bin1[covr.bin]
 binmed1<-binmed[covr.bin]}

covr.cat<-rep(F,length(catmed))
if(!is.null(catmed))
{for (i in 1:length(catmed))
 {if(testtype==1)
   covr.cat[i]<-ifelse(type3[xnames3==xname[catmed[i]],3]<alpha,T,F) #version 6: change to fit names instead of line numbers
  else if(testtype==2)
  {temp.data<-x[,c(catmed[i],pred)]
   if(surv)
    temp.fullmodel1<-coxph(y~.,data=data.frame(temp.data))
   else
    temp.fullmodel1<-glm(y~.,data=data.frame(temp.data),family=family1) #use type three error to test the full model
  temp.type3<-Anova(temp.fullmodel1,type="III")
  temp.p1<-temp.type3[rownames(temp.type3)==xname[catmed[i]],3]
  P1[grep(xname[catmed[i]],xnames3),]<-temp.p1
  covr.cat[i]<-ifelse(temp.p1<alpha,T,F)
  }
 }
  #covr.cat[i]<-ifelse(min(fullmodel$coef[grep(xname[catmed[i]],rownames(fullmodel$coef)),4])<alpha,T,F) #version 6: change to fit names instead of line numbers
 covr.cat<-ifelse(covr.cat+cat1>0,T,F)
 cat2<-cat1[covr.cat]
 catmed1<-catmed[covr.cat]
 catref1<-catref[covr.cat]}

a1<-c(contmed,binmed,catmed)
a2<-c(covr.cont,covr.bin,covr.cat)

cutx<-a1[!a2]

if (sum(a2)==0)
  return ("no mediators found")
else if(length(cutx)==0)
{newx1<-x
 contm1<-contmed
 binm1<-binmed
 catm1<-catmed
 catref1<-catref
 pred1<-pred}
else {newx1<-x[,-cutx]
      if(sum(covr.cont)==0)
      {contm1<-NULL}  
      else 
      {contm1<-colnum(contmed[covr.cont],cutx)}
      if(sum(covr.bin)==0)
        binm1<-NULL
      else 
        binm1<-colnum(binmed[covr.bin],cutx)
      if(sum(covr.cat)==0)
      {catm1<-NULL
       catref1<-NULL}
      else 
      {catm1<-colnum(catmed[covr.cat],cutx)
       catref1<-catref[covr.cat]}
      pred1<-colnum(pred,cutx)
      if(!is.null(jointm))
        for(i in 2:(jointm[[1]]+1))
          jointm[[i]]<-colnum(jointm[[i]],cutx)
}

#delete nonmediators
rela_var<-NULL               #to store the relationships between mediators and predictor   
rela_p<-NULL
name_newx<-names(newx1)
if (binpred)             #for binary (x)
{contm2<-contm1
 if(length(contm1)>0)
 {med.cont<-rep(F,length(contm1))
  for (i in 1:length(contm1))
  {tempmodel<-summary(glm(newx1[,contm1[i]]~newx1[,pred1]))
   med.cont[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }
  med.cont<-ifelse(med.cont+cont2>0,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
    {temp.p<-suppressWarnings(chisq.test(newx1[,pred1],newx1[,binm1[i]]))$p.value
     med.bin[i]<-ifelse(temp.p<alpha2,T,F)
     rela_var<-c(rela_var,name_newx[binm1[i]])
     rela_p<-c(rela_p,temp.p)
    }
  med.bin<-ifelse(med.bin+bin2>0,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))  
   {temp.p<-suppressWarnings(chisq.test(newx1[,pred1],newx1[,catm1[i]]))$p.value
    med.cat[i]<-ifelse(temp.p<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[catm1[i]])
    rela_p<-c(rela_p,temp.p)
   }
  med.cat<-ifelse(med.cat+cat2>0,T,F)
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}
else
{contm2<-contm1
 if(length(contm1)>0)
 {med.cont<-rep(F,length(contm1))
  for (i in 1:length(contm1))
  {temp.p<-cor.test(newx1[,contm1[i]],newx1[,pred1])$p.value
   med.cont[i]<-ifelse(temp.p<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-c(rela_p,temp.p)
  }
 med.cont<-ifelse(med.cont|cont2,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
  {tempmodel<-summary(glm(newx1[,pred1]~newx1[,binm1[i]]))
   med.bin[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[binm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }    
  med.bin<-ifelse(med.bin|bin2,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))   
  {tempmodel<-summary(glm(newx1[,pred1]~newx1[,catm1[i]]))
   med.cat[i]<-ifelse(tempmodel$coef[2,4]<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[catm1[i]])
   rela_p<-c(rela_p,tempmodel$coef[2,4])
  }
  med.cat<-(med.cat|cat2)
  cat3<-cat2[med.cat]
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}

if(length(catm2)==0)
  catm2<-NULL
if(length(binm2)==0)
  binm2<-NULL
if(length(contm2)==0)
  contm2<-NULL

newx2<-newx1

if (binpred) 
{if(!is.null(catm2))
  for (i in 1:length(catm2))
    newx2[,catm2[i]]<-as.factor(newx2[,catm2[i]])
 if(!is.null(binm2))
   for (i in 1:length(binm2))
     newx2[,binm2[i]]<-as.factor(newx2[,binm2[i]])
 results<-list(x=newx2, dirx=pred1, contm=contm2, catm=c(binm2,catm2),jointm=jointm,
              y=y,fullmodel=fullmodel1,rela=data.frame(name=rela_var,p=rela_p),binpred=binpred,surv=surv,
              testtype=testtype,P1=P1)
 class(results)<-"med_iden"
 return(results)
}
else
{catm<-NULL
 if(!is.null(catm2))
 {tempx<-cattobin(newx1,catm2,catref2)
  newx2<-tempx$x
  catm<-tempx$catm
  if(!is.null(binm2))
    binm2<-colnum(binm2,catm2)
  if(!is.null(contm2))
    contm2<-colnum(contm2,catm2)
  pred1<-colnum(pred1,catm2)
  if (!is.null(jointm))
  {if (sum(cat3)==0)
    for(i in 2:(jointm[[1]]+1))
      jointm[[i]]<-colnum(jointm[[i]],catm2)
   else
     for(i in 2:(jointm[[1]]+1))
     {a<-jointm[[i]]
      b<-NULL
      for (j in a)
        if(sum(j==catm2)==0)
          b<-c(b,colnum(j,catm2))
      else for (k in 1:length(catm2))
        if (j==catm2[k])
          b<-c(b,catm[[k+1]])
      jointm[[i]]<-b}
  }
 }
 results<-list(x=newx2,dirx=pred1,contm=contm2,binm=binm2,catm=catm, jointm=jointm, y=y,
               fullmodel=fullmodel1,rela=data.frame(name=rela_var,p=rela_p),binpred=binpred,surv=surv,
               testtype=testtype,P1=P1)
 class(results)<-"med_iden"
 return(results)
}
}

summary.med_iden<-function(object,...)
{var.name<-names(object$x)
 if(is.list(object$catm))                  #revised to show catm when it is a list (when x is continuous)
 {t.catm<-NULL
  for (i in 2:length(object$catm))
  t.catm<-c(t.catm,object$catm[[i]])}
 else
 {t.catm<-object$catm}
 mediator<-var.name[c(object$contm,object$binm,t.catm)]
 covariates<-var.name[-c(object$dirx,object$contm,object$binm,t.catm)]
 if(object$testtype==1)
  {tests<-Anova(object$fullmodel,type="III")
   temp<-rownames(tests)
   tests<-cbind(tests[,3],rep(NA,nrow(tests)))   
  }
 else
  {tests<-object$P1 
   temp<-rownames(tests)
   tests<-cbind(tests,rep(NA,nrow(tests)))}
  rownames(tests)<-temp
 for (i in 1:nrow(object$rela))
   tests[temp==object$rela[i,"name"],2]<-object$rela[i,"p"]
 dimnames(tests)[[2]]<-c("P-Value 1", "P-Value 2")
 result<-list(mediator=mediator, covariates=covariates,tests=tests, results=object)
 class(result)<-"summary.med_iden"
 result
}

print.summary.med_iden<-function(x,...)  #version 6: changed typos in the function--tests->x$tests
{cat("Identified as mediators: \n")
 print(x$mediator)
 cat("Selected as covariates: \n")
 print(x$covariates)
 cat("Tests: \n")
 temp<-rownames(x$tests)
 for (z in 1:length(temp))
   if(length(grep(temp[z],x$mediator))>0)
     dimnames(x$tests)[[1]][z]<-paste(temp[z],"*")
 if (!is.null(x$results$jointm))
  {tt<-NULL
   for (i in 1:(length(x$results$jointm)-1)) 
     tt<-c(tt,x$results$jointm[[i+1]])
   tt<-unique(tt)
   for (z in 1:length(temp))
     if(length(grep(temp[z],names(x$results$x)[tt]))>0)
       dimnames(x$tests)[[1]][z]<-paste(temp[z],"-")}
   #for (z in names(x$results$x)[tt])
     #rownames(x$tests)[grep(z,temp)]<-paste(temp[grep(z,temp)],"-")
 dimnames(x$tests)[[2]]<-c("P-Value 1", "P-Value 2")
 print(round(x$tests,3))
 cat("----\n *:mediator,-:joint mediator\n P-Value 1:Type-3 tests in the full model\n P-Value 2:Tests of relationship with the Predictor")
}

med<-function(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm,contm = data$contm, 
              catm = data$catm, jointm = data$jointm, allm = c(contm, catm), margin=1,
              n=20,seed=sample(1:1000,1), nonlinear=F, df=1, nu=0.001,D=3,
              distn=NULL,family1=NULL,refy=0,binpred=data$binpred,x.new=x,type=NULL)
{#for binary predictor
 med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                    catm = data$catm, jointm = data$jointm, 
                    allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                    D=3,distn="bernoulli",family1=binomial("logit"),
                    biny=F,refy=0,surv=F,type=NULL)
{if (is.null(allm))
  stop("Error: no potential mediator is specified")
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  
  if(biny)                     #recode y if y is binary
    y<-ifelse(y==refy,0,1)
  
  x<-x[!is.na(y),]             #delete nas in y for mart
  y<-y[!is.na(y)]
  
  te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
  {if(surv & !is.null(best.iter1))
    te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
  else if (surv)
    te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
  else
    te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
  te
  }
  
  med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
  {n3<-dim(x)[1]
  marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
  marg.m<-sample(marg.m)
  new1<-nom1
  new1[,med]<-marg.m[1:floor(n3/2)]
  new0<-nom0
  new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
  if(surv & !is.null(best.iter1))
    dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
  else if(surv)
    dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
  else
    dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
  dir.nom
  }
  
  med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
  {n3<-dim(x)[1]
  if (length(med)==1)                       #added for the new program, when there is only one mediator
   {if(is.factor(nom1[,med]))              #added to control for one factor mediator
     marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                as.character(sample(nom0[,med],replace=T)))))
    else
     marg.m<-sample(c(sample(nom1[,med],replace=T),
                    sample(nom0[,med],replace=T)))
   marg.m<-sample(marg.m)}        #added for the new program
  else                                         #added for the new program
   {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
  new1<-nom1
  new0<-nom0
  if(length(med)==1)                                       #added for the new program, when there is only one mediator
  {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
  new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
  else                                                     #added for the new program
  {new1[,med]<-marg.m[1:floor(n3/2),]
  new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
  if(surv & !is.null(best.iter1))
    dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
  else if(surv)
    dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
  else
    dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
  dir.nom
  }
  
  med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
  {n3<-dim(x)[1]
  marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
  marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
  dir.nom<-0
  for (i in levels(x[,med]))
  {new1<-nom1
  new1[1:dim(new1)[1],med]<-i
  new0<-nom0
  new0[1:dim(new0)[1],med]<-i
  p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
  if(surv & !is.null(best.iter1))
    dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
  else if(surv)
    dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
  else
    dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
  }
  dir.nom
  }
  
  #1.fit the model
  if (nonlinear)
  {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                        distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
  best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
  while(full.model$n.trees-best.iter1<30){
    full.model<-gbm.more(full.model, 50)           # do another 50 iterations
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
  else
  {if(surv)
    full.model<-coxph(y~., data=x)
   else
    full.model<-glm(y~., data=x, family=family1)
  best.iter1=NULL}
  #2. prepare for the store of results
  set.seed(seed)
  te<-rep(0,n)
  if(!is.null(jointm))
  {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
  dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  else
  {denm<-matrix(0,n,1+length(c(contm,catm)))
  dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
  ie<-denm
  #3. repeat to get the mediation effect
  for (k in 1:n)
  {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
    x1<-x[x[,dirx]==1,]
    x0<-x[x[,dirx]==0,]
    n3<-dim(x)[1]
    new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
    new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
    te[k]<-te.binx(full.model,new1,new0,best.iter1,surv,type)
    denm[k,1]<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
    j<-2
    #3.2 mediation effect from the continuous mediator
    if (!is.null(contm))
      for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type)
      j<-j+1}
    #3.3.mediation effect from the categorical mediator
    if (!is.null(catm))
      for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type)
      j<-j+1}
    #3.4 mediation effect from the joint mediators
    if (!is.null(jointm))
      for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[k,j]<-med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type)
      j<-j+1}
    
    #3.5 get the indirect effects
    ie[k,]<-te[k]-denm[k,]
    if(!is.null(jointm))
      dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
    else
      dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
  }
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
  class(a)<-"med"
  return(a)
}

#for continous predictor
 med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                    catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                    nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                    biny=F,refy=0,x.new=x,surv=F,type=NULL)
{if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(binm))
    binm<-unlist(sapply(binm,grep,xnames))
  if(!is.null(catm))
    for (i in 2:length(catm))
      if(is.character(catm[[i]]))
        catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  if(biny)                     #recode y if y is binary
    y<-ifelse(y==refy,0,1)
  
  x<-x[!is.na(y),]             #delete nas in y for mart
  y<-y[!is.na(y)]
  
  anymissing<-function(vec)
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
  col_mean<-function(col,n.row)
  {temp<-matrix(col,n.row)
  apply(temp,1,mean,na.rm=T)}
  
  
  dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
  {models<-NULL
  res<-NULL
  if(!is.null(catm))
  {for (i in 2:(catm$n+1))
    binm<-c(binm,catm[[i]])}
  
  z<-x[,dirx]
  j<-1
  if(!is.null(binm))
  {for(i in binm)
  {if(nonlinear)
    models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
  else
    models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
  res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
  j<-j+1}
  }
  for (i in contm)
  {if(nonlinear)
    models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
  else
    models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
  res<-cbind(res,models[[j]]$res)
  j<-j+1
  }
  list(models=models,varmat=var(res))
  }
  
  
  sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
  {mult.norm<-function(mu,vari,n) 
  {  if (nrow(vari)!=ncol(vari)) 
  {result<-c("Error: Variance matrix is not square")}  
    else if (length(mu)!=nrow(vari)) 
    {result<-c("Error: length mu is not right!")}  
    else {   p<-length(mu)
    tmp1<-eigen(vari)$values
    tmp2<-eigen(vari)$vectors   
    result<-matrix(0,n,p)   
    for (i in 1:p)
    {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
    for (i in 1:n)
    {result[i,]<-tmp2%*%result[i,]+mu}
    }  
    result
  }
  
  match.margin<-function(vec)   
  {range1<-vec[1:2]
  vec1<-vec[-(1:2)]
  range2<-range(vec1,na.rm=T)
  vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
  vec1
  }
  
  gen.mult<-function(vec)
  {if(sum(is.na(vec))>0)
    return(rep(NA,length(vec)))
    else{ 
      l<-1-sum(vec)
      l<-ifelse(l<0,0,l)
      return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
  }
  
  means<-NULL
  z<-x1[,dirx]
  binm1<-binm
  if(!is.null(catm))
  {for (i in 2:(catm$n+1))
    binm1<-c(binm1,catm[[i]])}
  
  if(!is.null(binm1))
    for (i in 1:length(binm1))
      means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
  
  if(!is.null(contm))
    for (i in (length(binm1)+1):length(c(binm1,contm)))
      means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
  if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
  {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
  sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
  else{
    sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
    
    range.means<-apply(means,2,range,na.rm=T)
    
    sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
  }
  sim.m2<-data.frame(sim.m2)
  n<-dim(sim.m2)[1]
  if(!is.null(binm))
    for (i in 1:length(binm))
      sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
  
  if(!is.null(catm))
  {j<-length(binm)+1
  for (i in 2:(catm$n+1))
  {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
  sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
  j<-j+length(catm[[i]])}
  }
  
  x1[,c(binm1,contm)]<-sim.m2
  
  x1
  }
  
  
  if(is.null(catm))
    multi=jointm
  else if(is.null(jointm))
    multi=catm
  else {temp1<-catm
  temp2<-jointm
  temp1[[1]]=catm[[1]]+jointm[[1]]
  temp2[[1]]<-NULL
  multi=append(temp1,temp2)} 
  listm=list(single=c(contm,binm),multi=multi)
  
  nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
  x<-x[nonmissing,]
  y<-y[nonmissing]
  nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
  x.new<-x.new[nonmissing1,]
  
  #1.fit the model
  if(nonlinear)
  {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                        distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
  best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
  while(full.model$n.trees-best.iter1<30){
    full.model<-gbm.more(full.model, 50)           # do another 50 iterations
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
  }
  else
  {if(surv)
    full.model<-coxph(y~., data=x)
   else
    full.model<-glm(y~., data=x, family=family1)
  best.iter1=NULL}
  
  #2. prepare for the store of results
  set.seed(seed)
  n.new<-nrow(x.new)
  te<-rep(0,n.new)
  
  #3. get the joint distribution of m given x
  distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
  te1<-NULL
  #x1<-x.new
  #x1[,dirx]<-x[,dirx]+margin
  #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
  
  n1<-dim(x)[1]
  denm1<-NULL

  #4. repeat to get the mediation effect
  for (k in 1:n)
  {new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
  x1<-new0
  x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
  if(surv & !is.null(best.iter1))
  {ybar0<-predict(full.model,new0,best.iter1,type=type)
   denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
  else if(surv)
  {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
   denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
  else
    {ybar0<-predict(full.model,new0,best.iter1)
     denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
  new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
  #4.1 get the te
  if(surv & !is.null(best.iter1))
    te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type)-ybar0)/margin)
  else if(surv)
    te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
  else
    te1<-cbind(te1, (predict(full.model,new1,best.iter1)-ybar0)/margin)
  
  #4.2 mediation effect from the single mediator
  if (!is.null(listm$single))
    for (i in 1:length(listm$single))
    {new1.nm<-new1
    new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
    if(surv & !is.null(best.iter1))
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
    }
  
  #4.3.mediation effect from the joint mediator
  if (!is.null(listm$multi))
    for (i in 2:(listm$multi[[1]]+1))
    {new1.nm<-new1
    new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
    if(surv & !is.null(best.iter1))
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
    }
  denm1<-rbind(denm1,denm2)
  } 
  
  #4.4 get the indirect effects
  denm<-apply(denm1,2,col_mean,n.new)
  te<-apply(te1,1,mean)
  ie<-te-denm
  if(!is.null(listm$multi))
    colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
  else 
    colnames(denm)<-c("de",names(x)[listm$single])
  if(!is.null(listm$multi))
    colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
  else 
    colnames(ie)<-c("all",names(x)[listm$single])
  
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new)
  class(a)<-"med"
  return(a)
}

 surv=F
 if(class(y)=="Surv")
 {surv=T
  biny=F
  if(is.null(distn))
   distn = "coxph"
  if(is.null(type) & nonlinear)
    type="response"
  else if (is.null(type))
    type="risk"}
 else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
  {biny=T
   if(is.null(family1))
     family1 = binomial("logit")
   if(is.null(distn))
     distn = "bernoulli"}
  else
  {biny=F
   if(is.null(family1))
    family1 = gaussian(link = "identity")
   if(is.null(distn))
    distn = "gaussian"}
  
  if(binpred)
    a<-med.binx(data=data, x=x, y=y, dirx=dirx, contm = contm, 
                catm=catm, jointm=jointm, allm=allm, n=n,seed=seed,
                nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,
                biny=biny,refy=refy,surv=surv,type=type)
  else
    a<-med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                 catm=catm, jointm=jointm, margin=margin, n=n, seed=seed, 
                 nonlinear=nonlinear, df=df, nu=nu,D=D, distn=distn, 
                 family1=family1,biny=biny,refy=refy,x.new=x.new,surv=surv,type=type)
  return(a)
}
  

print.med<-function(x,...,digit=4)
{cat("The estimated total effect:")
  print(round(mean(x$te,na.rm=T),digit))
 cat("\nThe estimated indirect effect:\n")
  print(round(apply(x$ie,2,mean,na.rm=T),digit))
}

boot.med<-function(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,catm=data$catm,
                        jointm=data$jointm,margin=1,n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,
                        D=3,distn=NULL,family1=NULL,n2=50,
                        weight=rep(1,nrow(x)),refy=NULL,x.new=x,binpred=data$binpred,type=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         weight=rep(1,nrow(x)),biny=F,refy=0,surv,type)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, 
                     allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn="bernoulli",family1=binomial("logit"),
                     biny=F,refy=0,surv=F,type=NULL)
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {if(surv & !is.null(best.iter1))
      te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if (surv)
      te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    te
    }
    
    med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
    marg.m<-sample(marg.m)
    new1<-nom1
    new1[,med]<-marg.m[1:floor(n3/2)]
    new0<-nom0
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                 as.character(sample(nom0[,med],replace=T)))))
    else
      marg.m<-sample(c(sample(nom1[,med],replace=T),
                       sample(nom0[,med],replace=T)))
    marg.m<-sample(marg.m)}        #added for the new program
    else                                         #added for the new program
    {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
    new1<-nom1
    new0<-nom0
    if(length(med)==1)                                       #added for the new program, when there is only one mediator
    {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
    else                                                     #added for the new program
    {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
    marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
    dir.nom<-0
    for (i in levels(x[,med]))
    {new1<-nom1
    new1[1:dim(new1)[1],med]<-i
    new0<-nom0
    new0[1:dim(new0)[1],med]<-i
    p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
    if(surv & !is.null(best.iter1))
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
    else if(surv)
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
    else
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
    }
    dir.nom
    }
    
    #1.fit the model
    if (nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    #2. prepare for the store of results
    set.seed(seed)
    te<-rep(0,n)
    if(!is.null(jointm))
    {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
    else
    {denm<-matrix(0,n,1+length(c(contm,catm)))
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {    x1<-x[x[,dirx]==1,]
    x0<-x[x[,dirx]==0,]
    n3<-dim(x)[1]
    new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
    new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
    #3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      te[k]<-te.binx(full.model,new1,new0,best.iter1,surv,type)
      denm[k,1]<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type)
        j<-j+1}
      
      #3.5 get the indirect effects
      ie[k,]<-te[k]-denm[k,]
      if(!is.null(jointm))
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
    }
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  if(is.null(jointm))
  {ie<-matrix(0,n2+1,1+length(c(contm,catm)))
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])}
  else
  {ie<-matrix(0,n2+1,1+length(c(contm,catm))+jointm[[1]])
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  temp<-med.binx(data=F,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1]<-mean(temp$te)
  de[1]<-mean(temp$denm[,1])
  ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  model<-temp$model
  best.iter<-temp$best.iter
  for (i in 1:n2)
  {boots<-sample(1:nrow(x),replace=T,prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.binx(data=F,x1,y1,dirx,contm,catm,jointm,allm,n,seed+i,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1+i]<-mean(temp$te)
  de[1+i]<-mean(temp$denm[,1])
  ie[1+i,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  print(i)
  }
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         weight=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                      biny=F,refy=0,x.new=x,surv=F,type=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    anymissing<-function(vec)
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row)
    {temp<-matrix(col,n.row)
    apply(temp,1,mean,na.rm=T)}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-x[,dirx]
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
    else
      models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
    {mult.norm<-function(mu,vari,n) 
    {  if (nrow(vari)!=ncol(vari)) 
    {result<-c("Error: Variance matrix is not square")}  
      else if (length(mu)!=nrow(vari)) 
      {result<-c("Error: length mu is not right!")}  
      else {   p<-length(mu)
      tmp1<-eigen(vari)$values
      tmp2<-eigen(vari)$vectors   
      result<-matrix(0,n,p)   
      for (i in 1:p)
      {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
      for (i in 1:n)
      {result[i,]<-tmp2%*%result[i,]+mu}
      }  
      result
    }
    
    match.margin<-function(vec)   
    {range1<-vec[1:2]
    vec1<-vec[-(1:2)]
    range2<-range(vec1,na.rm=T)
    vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
    vec1
    }
    
    gen.mult<-function(vec)
    {if(sum(is.na(vec))>0)
      return(rep(NA,length(vec)))
      else{ 
        l<-1-sum(vec)
        l<-ifelse(l<0,0,l)
        return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
    }
    
    means<-NULL
    z<-x1[,dirx]
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
    if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
    {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
    sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
    else{
      sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
      
      range.means<-apply(means,2,range,na.rm=T)
      
      sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
    }
    sim.m2<-data.frame(sim.m2)
    n<-dim(sim.m2)[1]
    if(!is.null(binm))
      for (i in 1:length(binm))
        sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
    
    if(!is.null(catm))
    {j<-length(binm)+1
    for (i in 2:(catm$n+1))
    {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
    sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
    j<-j+length(catm[[i]])}
    }
    
    x1[,c(binm1,contm)]<-sim.m2
    
    x1
    }
    
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
    x<-x[nonmissing,]
    y<-y[nonmissing]
    nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
    x.new<-x.new[nonmissing1,]
    
    #1.fit the model
    if(nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
    }
    else
    {if(surv)
      full.model<-suppressWarnings(coxph(y~., data=x))
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    te<-rep(0,n.new)
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
    te1<-NULL
    #x1<-x.new
    #x1[,dirx]<-x[,dirx]+margin
    #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
    
    n1<-dim(x)[1]
    denm1<-NULL
    
    #4. repeat to get the mediation effect
    for (k in 1:n)
    {new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
    x1<-new0
    x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
    if(surv & !is.null(best.iter1))
    {ybar0<-predict(full.model,new0,best.iter1,type=type)
    denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
    else if(surv)
    {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
    denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
    else
    {ybar0<-predict(full.model,new0,best.iter1)
    denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
    new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
    #4.1 get the te
    if(surv & !is.null(best.iter1))
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      te1<-cbind(te1, (predict(full.model,new1,best.iter1)-ybar0)/margin)
    
    #4.2 mediation effect from the single mediator
    if (!is.null(listm$single))
      for (i in 1:length(listm$single))
      {new1.nm<-new1
      new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
      for (i in 2:(listm$multi[[1]]+1))
      {new1.nm<-new1
      new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    denm1<-rbind(denm1,denm2)
    } 
    
    #4.4 get the indirect effects
    denm<-apply(denm1,2,col_mean,n.new)
    te<-apply(te1,1,mean)
    ie<-te-denm
    if(!is.null(listm$multi))
      colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(denm)<-c("de",names(x)[listm$single])
    if(!is.null(listm$multi))
      colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(ie)<-c("all",names(x)[listm$single])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new)
    class(a)<-"med"
    return(a)
  }
  
if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")

xnames<-colnames(x)
if(is.character(dirx))
  pred<-grep(dirx,xnames)
if(is.character(contm))
  contm<-unlist(sapply(contm,grep,xnames))
if(is.character(binm))
  binm<-unlist(sapply(binm,grep,xnames))
if(!is.null(catm))
  for (i in 2:length(catm))
    if(is.character(catm[[i]]))
      catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
if(!is.null(jointm))
  for (i in 2:length(jointm))
    if(is.character(jointm[[i]]))
      jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))

if(is.null(catm))
{multi=jointm
name1<-NULL                       #added in the new program
if (!is.null(multi))              #added in the new program, in case that multi is NULL
  name1<-paste("j",1:multi[[1]],sep="")}
else if(is.null(jointm))
{multi=catm
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,names(x)[multi[[i]][1]])}
else {temp1<-catm
temp2<-jointm
temp1[[1]]=catm[[1]]+jointm[[1]]
temp2[[1]]<-NULL
multi=append(temp1,temp2)
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,names(x)[multi[[i]][1]])
name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
listm=list(single=c(contm,binm),multi=multi)

te<-rep(0,n2+1)
de<-rep(0,n2+1)
mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
ie<-matrix(0,n2+1,1+length(listm$single)+mul)   #added in the new program
if(!is.null(listm$multi))
  dimnames(ie)[[2]]<-c("all",names(x)[listm$single],name1)
else 
  dimnames(ie)[[2]]<-c("all",names(x)[listm$single])

temp<-med.contx(data=F,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                refy=refy,x.new=x.new,surv,type)
x.new<-temp$x.new
te[1]<-mean(temp$te,na.rm=T)
de[1]<-mean(temp$denm[,1],na.rm=T) 
ie[1,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
te1<-NULL                      #to store the mediation effects on predictor
de1<-NULL
ie1<-NULL
model<-temp$model
for (l in 1:n2)
{boots<-sample(1:nrow(x),replace=T, prob=weight)
x1<-x[boots,]
y1<-y[boots]
temp<-med.contx(data=F,x1,y1,dirx,binm,contm,catm,jointm, margin,n,seed+l,nonlinear,df,nu,D,distn,family1,biny,refy,x.new,surv,type) #added to the new codel, change the seed to make different results
te[1+l]<-mean(temp$te,na.rm=T)
de[1+l]<-mean(temp$denm[,1],na.rm=T)
ie[1+l,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
te1<-cbind(te1,temp$te)
de1<-cbind(de1,temp$denm[,1])
ie1<-rbind(ie1,temp$ie)
print(l)
}
a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model,
        data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
        boot.detail=list(x.new=x.new[,dirx],te1=te1,de1=de1,ie1=ie1))
class(a)<-"mma"
return(a)
}

surv=F
if(class(y)=="Surv")
{surv=T
 biny=F
 if(is.null(distn))
  distn = "coxph"
 if(is.null(type) & nonlinear)
  type="response"
 else if (is.null(type))
  type="risk"}
else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
{biny=T
if(is.null(family1))
  family1 = binomial("logit")
if(is.null(distn))
  distn = "bernoulli"
if(!is.null(refy))
  y<-ifelse(y==refy,0,1)
else
  y<-ifelse(as.factor(y)==levels(as.factor(y))[1],0,1)
}
else
{biny=F
if(is.null(family1))
  family1 = gaussian(link = "identity")
if(is.null(distn))
  distn = "gaussian"
}

if(binpred)
  a<-boot.med.binx(data=data,x=x, y=y,dirx=dirx,contm=contm,catm=catm,
                             jointm=jointm,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,
                             D=D,distn=distn,family1=family1,
                             weight=weight,biny=biny,refy=0,surv,type)
else
  a<-boot.med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                    catm=catm, jointm=jointm, margin = margin, n = n, seed = seed, 
                    nonlinear = nonlinear, df = df, nu = nu, D = D, distn = distn, 
                    family1 = family1, n2 = n2,weight=weight,
                    biny=biny,refy=0,x.new=x.new,surv,type)

return(a)
}

    
  


mma<-function(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,
              catmed=NULL,catref=NULL,jointm=NULL,refy=NULL,
              predref=NULL,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
              nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=NULL,n2=50,weight=rep(1,nrow(x)),x.new=NULL,type=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         weight=rep(1,nrow(x)),biny=F,refy=0,surv,type)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, 
                     allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn="bernoulli",family1=binomial("logit"),
                     biny=F,refy=0,surv=F,type=NULL)
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {if(surv & !is.null(best.iter1))
      te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if (surv)
      te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    te
    }
    
    med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
    marg.m<-sample(marg.m)
    new1<-nom1
    new1[,med]<-marg.m[1:floor(n3/2)]
    new0<-nom0
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                 as.character(sample(nom0[,med],replace=T)))))
    else
      marg.m<-sample(c(sample(nom1[,med],replace=T),
                       sample(nom0[,med],replace=T)))
    marg.m<-sample(marg.m)}        #added for the new program
    else                                         #added for the new program
    {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
    new1<-nom1
    new0<-nom0
    if(length(med)==1)                                       #added for the new program, when there is only one mediator
    {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
    else                                                     #added for the new program
    {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
    marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
    dir.nom<-0
    for (i in levels(x[,med]))
    {new1<-nom1
    new1[1:dim(new1)[1],med]<-i
    new0<-nom0
    new0[1:dim(new0)[1],med]<-i
    p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
    if(surv & !is.null(best.iter1))
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
    else if(surv)
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
    else
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
    }
    dir.nom
    }
    
    #1.fit the model
    if (nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    #2. prepare for the store of results
    set.seed(seed)
    te<-rep(0,n)
    if(!is.null(jointm))
    {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
    else
    {denm<-matrix(0,n,1+length(c(contm,catm)))
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {    x1<-x[x[,dirx]==1,]
    x0<-x[x[,dirx]==0,]
    n3<-dim(x)[1]
    new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
    new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
    #3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      te[k]<-te.binx(full.model,new1,new0,best.iter1,surv,type)
      denm[k,1]<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type)
        j<-j+1}

            #3.5 get the indirect effects
      ie[k,]<-te[k]-denm[k,]
      if(!is.null(jointm))
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
    }
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  if(is.null(jointm))
  {ie<-matrix(0,n2+1,1+length(c(contm,catm)))
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])}
  else
  {ie<-matrix(0,n2+1,1+length(c(contm,catm))+jointm[[1]])
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  temp<-med.binx(data=F,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1]<-mean(temp$te)
  de[1]<-mean(temp$denm[,1])
  ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  model<-temp$model
  best.iter<-temp$best.iter
  for (i in 1:n2)
  {boots<-sample(1:nrow(x),replace=T,prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.binx(data=F,x1,y1,dirx,contm,catm,jointm,allm,n,seed+i,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1+i]<-mean(temp$te)
  de[1+i]<-mean(temp$denm[,1])
  ie[1+i,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  print(i)
  }
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         weight=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                      biny=F,refy=0,x.new=x,surv=F,type=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    anymissing<-function(vec)
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row)
    {temp<-matrix(col,n.row)
    apply(temp,1,mean,na.rm=T)}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-x[,dirx]
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
    else
      models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
    {mult.norm<-function(mu,vari,n) 
    {  if (nrow(vari)!=ncol(vari)) 
    {result<-c("Error: Variance matrix is not square")}  
      else if (length(mu)!=nrow(vari)) 
      {result<-c("Error: length mu is not right!")}  
      else {   p<-length(mu)
      tmp1<-eigen(vari)$values
      tmp2<-eigen(vari)$vectors   
      result<-matrix(0,n,p)   
      for (i in 1:p)
      {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
      for (i in 1:n)
      {result[i,]<-tmp2%*%result[i,]+mu}
      }  
      result
    }
    
    match.margin<-function(vec)   
    {range1<-vec[1:2]
    vec1<-vec[-(1:2)]
    range2<-range(vec1,na.rm=T)
    vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
    vec1
    }
    
    gen.mult<-function(vec)
    {if(sum(is.na(vec))>0)
      return(rep(NA,length(vec)))
      else{ 
        l<-1-sum(vec)
        l<-ifelse(l<0,0,l)
        return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
    }
    
    means<-NULL
    z<-x1[,dirx]
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
    if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
    {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
    sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
    else{
      sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
      
      range.means<-apply(means,2,range,na.rm=T)
      
      sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
    }
    sim.m2<-data.frame(sim.m2)
    n<-dim(sim.m2)[1]
    if(!is.null(binm))
      for (i in 1:length(binm))
        sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
    
    if(!is.null(catm))
    {j<-length(binm)+1
    for (i in 2:(catm$n+1))
    {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
    sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
    j<-j+length(catm[[i]])}
    }
    
    x1[,c(binm1,contm)]<-sim.m2
    
    x1
    }
    
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
    x<-x[nonmissing,]
    y<-y[nonmissing]
    nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
    x.new<-x.new[nonmissing1,]
   
    #1.fit the model
    if(nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
    }
    else
    {if(surv)
      full.model<-suppressWarnings(coxph(y~., data=x))
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    te<-rep(0,n.new)
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
    te1<-NULL
    #x1<-x.new
    #x1[,dirx]<-x[,dirx]+margin
    #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
    
    n1<-dim(x)[1]
    denm1<-NULL
    
    #4. repeat to get the mediation effect
    for (k in 1:n)
    {new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
    x1<-new0
    x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
    if(surv & !is.null(best.iter1))
    {ybar0<-predict(full.model,new0,best.iter1,type=type)
    denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
    else if(surv)
    {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
    denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
    else
    {ybar0<-predict(full.model,new0,best.iter1)
    denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
    new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
    #4.1 get the te
    if(surv & !is.null(best.iter1))
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      te1<-cbind(te1, (predict(full.model,new1,best.iter1)-ybar0)/margin)
    
    #4.2 mediation effect from the single mediator
    if (!is.null(listm$single))
      for (i in 1:length(listm$single))
      {new1.nm<-new1
      new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
      for (i in 2:(listm$multi[[1]]+1))
      {new1.nm<-new1
      new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    denm1<-rbind(denm1,denm2)
    } 
    
    #4.4 get the indirect effects
    denm<-apply(denm1,2,col_mean,n.new)
    te<-apply(te1,1,mean)
    ie<-te-denm
    if(!is.null(listm$multi))
      colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(denm)<-c("de",names(x)[listm$single])
    if(!is.null(listm$multi))
      colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(ie)<-c("all",names(x)[listm$single])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type,  model=full.model,best.iter=best.iter1),x.new=x.new)
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(binm))
    binm<-unlist(sapply(binm,grep,xnames))
  if(!is.null(catm))
    for (i in 2:length(catm))
      if(is.character(catm[[i]]))
        catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  if(is.null(catm))
  {multi=jointm
  name1<-NULL                       #added in the new program
  if (!is.null(multi))              #added in the new program, in case that multi is NULL
    name1<-paste("j",1:multi[[1]],sep="")}
  else if(is.null(jointm))
  {multi=catm
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])}
  else {temp1<-catm
  temp2<-jointm
  temp1[[1]]=catm[[1]]+jointm[[1]]
  temp2[[1]]<-NULL
  multi=append(temp1,temp2)
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])
  name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
  listm=list(single=c(contm,binm),multi=multi)
  
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2+1,1+length(listm$single)+mul)   #added in the new program
  if(!is.null(listm$multi))
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single],name1)
  else 
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single])
  
  temp<-med.contx(data=F,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,surv,type)
  x.new=temp$x.new
  te[1]<-mean(temp$te,na.rm=T)
  de[1]<-mean(temp$denm[,1],na.rm=T) 
  ie[1,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  te1<-NULL                      #to store the mediation effects on predictor
  de1<-NULL
  ie1<-NULL
  model<-temp$model
  for (l in 1:n2)
  {boots<-sample(1:nrow(x),replace=T, prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.contx(data=F,x1,y1,dirx,binm,contm,catm,jointm, margin,n,seed+l,nonlinear,df,nu,D,distn,family1,biny,refy,x.new,surv,type) #added to the new codel, change the seed to make different results
  te[1+l]<-mean(temp$te,na.rm=T)
  de[1+l]<-mean(temp$denm[,1],na.rm=T)
  ie[1+l,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  te1<-cbind(te1,temp$te)
  de1<-cbind(de1,temp$denm[,1])
  ie1<-rbind(ie1,temp$ie)
  print(l)
  }
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(x.new=x.new[,dirx],te1=te1,de1=de1,ie1=ie1))
  class(a)<-"mma"
  return(a)
}


surv=F
if(class(y)=="Surv")
{surv=T
biny=F
if(is.null(distn))
  distn = "coxph"
if(is.null(type) & nonlinear)
  type="response"
else if (is.null(type))
  type="risk"}
else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
{biny=T
if(is.null(family1))
  family1 = binomial("logit")
if(is.null(distn))
  distn = "bernoulli"
if(!is.null(refy))
  y<-ifelse(y==refy,0,1)
else
  y<-ifelse(as.factor(y)==levels(as.factor(y))[1],0,1)
}
else
{biny=F
if(is.null(family1))
  family1 = gaussian(link = "identity")
if(is.null(distn))
  distn = "gaussian"
}

data<-data.org(x=x,y=y,pred=pred,mediator=mediator,contmed=contmed,binmed=binmed,
               binref=binref,catmed=catmed,catref=catref,jointm=jointm,refy=refy,family1=family1,
               predref=predref,alpha=alpha,alpha2=alpha2)
binpred<-data$binpred

if(binpred) 
  {result<-boot.med.binx(data=data,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,weight=weight,biny=biny,refy=0,surv=surv,type=type)
  }
 else
  {if(is.null(x.new))
    result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                           D=D,distn=distn,family1=family1,n2=n2,weight=weight,biny=biny,refy=0,surv=surv,type=type)
   else
    result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                           D=D,distn=distn,family1=family1,n2=n2,weight=weight,biny=biny,refy=0, x.new=x.new,surv=surv,type=type)
   
  }
 result
}


#classes and methods for mma
print.mma<-function(x,...)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$model$MART)
   cat ("MART\n")
 else cat("GLM\n")
 print(x$e)
}


summary.mma<-function(object,...,alpha=0.05,plot=TRUE,RE=FALSE,quant=T)
{x<-object
 temp1<-x$boots
 temp2<-x$est
 temp3<-x$boots   #calculate the RE
 temp3$ie<-temp3$ie/temp3$te
 temp3$de<-temp3$de/temp3$te
 temp4<-x$est
 temp4$ie<-temp4$ie/temp4$te
 temp4$de<-temp4$de/temp4$te
 
 a1<-alpha/2
 a2<-1-a1
 b1<-qnorm(a1)
 b2<-qnorm(a2)
 temp1.result<-list(indirect.effect=rbind(est=temp2$ie,mean=apply(temp1$ie,2,mean),sd=apply(temp1$ie,2,sd),
                                          upbd=apply(temp1$ie,2,mean)+b2*apply(temp1$ie,2,sd),
                                          lwbd=apply(temp1$ie,2,mean)+b1*apply(temp1$ie,2,sd),
                                          upbd_q=apply(temp1$ie,2,quantile,a2), lwbd_q=apply(temp1$ie,2,quantile,a1)),
                    total.effect=c(est=temp2$te,mean=mean(temp1$te),sd=sd(temp1$te),
                                   upbd=mean(temp1$te)+b2*sd(temp1$te),lwbd=mean(temp1$te)+b1*sd(temp1$te),
                                   upbd=quantile(temp1$te,a2,na.rm=T),lwbd=quantile(temp1$te,a1,na.rm=T)),
                    direct.effect=c(est=temp2$de,mean=mean(temp1$de),sd=sd(temp1$de),
                                    upbd=mean(temp1$de)+b2*sd(temp1$de),lwbd=mean(temp1$de)+b1*sd(temp1$de),
                                    upbd=quantile(temp1$de,a2,na.rm=T),lwbd=quantile(temp1$de,a1,na.rm=T)))
 temp2.result<-list(indirect.effect=rbind(est=temp4$ie,mean=apply(temp3$ie,2,mean),sd=apply(temp3$ie,2,sd),
                                          upbd=apply(temp3$ie,2,mean)+b2*apply(temp3$ie,2,sd),
                                          lwbd=apply(temp3$ie,2,mean)+b1*apply(temp3$ie,2,sd),
                                          upbd_q=apply(temp3$ie,2,quantile,a2,na.rm=T), 
                                          lwbd_q=apply(temp3$ie,2,quantile,a1,na.rm=T)),
                    direct.effect=c(est=temp4$de,mean=mean(temp3$de),sd=sd(temp3$de),
                                    upbd=mean(temp3$de)+b2*sd(temp3$de),lwbd=mean(temp3$de)+b1*sd(temp3$de),
                                    upbd=quantile(temp3$de,a2,na.rm=T),lwbd=quantile(temp3$de,a1,na.rm=T)))
 result<-list(results=temp1.result,re=temp2.result,alpha=alpha,plot=plot,obj=x,RE=RE,quant=quant)
 class(result)<-"summary.mma"
 result
 }

print.summary.mma<-function(x,...)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$obj$model$MART)
  cat ("MART\n")
 else cat("GLM\n")
 print(lapply(x$results,round,3))
 if(x$RE)
 {cat("The relative effects:\n")
  print(lapply(x$re,round,3))}
if(x$plot)
{re<-c(x$re$indirect.effect[2,-1],x$re$dir[2])
 if(x$quant)
 {upper<-c(x$re$indirect.effect[6,-1],x$re$dir[6])
  lower<-c(x$re$indirect.effect[7,-1],x$re$dir[7])}
 else
 {upper<-c(x$re$indirect.effect[4,-1],x$re$dir[4])
  lower<-c(x$re$indirect.effect[5,-1],x$re$dir[5])}
 d<-order(re)
 name1<-c(colnames(x$re$indirect.effect)[-1],"de")
 par(mfrow=c(1,1),mar=c(1,6,1,1),oma=c(3,2,2,4))
 bp <- barplot2(re[d], horiz = TRUE, main="Relative Effects", 
                names.arg=name1[d],plot.ci = TRUE, ci.u = upper[d], ci.l = lower[d],
                cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower)),
                col = rainbow(length(d), start = 3/6, end = 4/6))
}
}


plot.mma<-function(x,...,vari,xlim=range(x$data$x[,vari],na.rm=T),alpha=0.95,quantile=F)
{marg.den<-function(x,y)
{y<-y[!is.na(x)]
 x<-x[!is.na(x)]
 z1<-unique(x)
 z2<-rep(0,length(z1))
 for (i in 1:length(z1))
   z2[i]<-mean(y[x==z1[i]],na.rm=T)
 z3<-order(z1)
 cbind(z1[z3],z2[z3])
}

overlapHist <- function(a, b,breaks=NULL, xlim=NULL, xname=NULL)
{a1<-a
b1<-b
a<-a[!is.na(a1) & !is.na(b1)]
b<-b[!is.na(a1) & !is.na(b1)]
j<-unique(b)
ahist<-hist(a[b==j[1]],plot=F)
dist = ahist$breaks[2]-ahist$breaks[1]
lb =min(ahist$breaks,na.rm = T)
ub=max(ahist$breaks,na.rm = T)
yl=max(ahist$density,na.rm = T)
for(i in j[-1])
{bhist<-hist(a[b==i],plot=F)
lb =min(lb,bhist$breaks,na.rm = T)
ub =max(ub,bhist$breaks,na.rm = T)
yl=max(yl,bhist$density,na.rm = T)
dist = min(dist,bhist$breaks[2]-bhist$breaks[1])
}
breaks=seq(lb,ub,dist)
if(is.null(xlim))
  xlim=c(lb,ub)
for (i in j)
  hist(a[b==i],ylab="Density",xlab="",breaks=breaks, 
       xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="="))
}


boot.ci<-function(x,mat,alpha,quantile=F) #the mat is the booted results with row be different x, and columns diff boot
  #cri_val is the critical value
{x.uniq<-sort(unique(x,na.rm=T))
mn<-NULL
upbd<-NULL
lwbd<-NULL
alpha<-(1-alpha)/2
for (i in x.uniq)
{sd_dev<-sd(as.vector(mat[x==i,]),na.rm=T)
mn1<-mean(as.vector(mat[x==i,]),na.rm=T)
if(quantile)
{upbd<-c(upbd,quantile(as.vector(mat[x==i,]),1-alpha,na.rm=T))
lwbd<-c(lwbd,quantile(as.vector(mat[x==i,]),alpha,na.rm=T))
}
else
{cri_val<-qnorm(1-alpha)
upbd<-c(upbd,mn1+cri_val*sd_dev)
lwbd<-c(lwbd,mn1-cri_val*sd_dev)}
mn<-c(mn,mn1)}
x.uniq<-x.uniq[!is.na(lwbd)&!is.na(upbd)]
tt<-(!is.na(lwbd)) & (!is.na(upbd))
mn<-mn[tt]
lwbd<-lwbd[tt]
upbd<-upbd[tt]
return(data.frame(x=x.uniq,F=mn,L=lwbd,U=upbd))
}
plot_ci<-function(df,xlab="x",ylab="IE")
{plot(df$x, df$F, ylim = range(c(df$L,df$U)), type = "l",xlab=xlab,ylab=ylab)
  polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
  lines(df$x, df$F, lwd = 2)
  lines(df$x, df$U, col="red",lty=2)
  lines(df$x, df$L, col="red",lty=2)}
op <- par(no.readonly = TRUE) # the whole list of settable par's.
 if (x$model[1]==T) 
 {full.model=x$model$model
  best.iter=x$model$best.iter
  data=x$data
  mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
  if(data$binpred)
  {if(!is.factor(data$x[,vari]))
   {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    overlapHist(a=data$x[,vari],b=data$x[,data$dirx],xlim=xlim,xname=names(data$x)[data$dirx])
    #ggplot(data$x[,vari], aes(length, fill = data$x[,data$dirx])) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    #for (j in unique(data$x[,data$dirx]))
    # hist(data$x[data$x[,data$dirx]==j,vari],xlab="",ylab="Relative Freq",xlim=xlim, freq=F,main=paste(names(data$x)[data$dirx],j,sep="="))
    }
   else{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
     if(full.model$distribution=="gaussian")
       suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
     else if(full.model$distribution=="coxph")
       suppressWarnings(plot.gbm(full.model, i.var=vari))
     else
          suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
     #plot(data$x[,vari]~data$x[,data$dirx],ylab=mname,xlab=names(data$x)[data$dirx])}
     temp1<-NULL
     for (j in unique(data$x[,data$dirx]))
       temp1<-c(temp1,prop.table(table(data$x[data$x[,data$dirx]==j,vari])))
     for (j in unique(data$x[,data$dirx]))
       barplot(prop.table(table(data$x[data$x[,data$dirx]==j,vari])),ylim=c(0,max(temp1,na.rm=T)),
               ylab="Prop",sub=paste(names(data$x)[data$dirx],j,sep="="))}
  }
  else
  {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    ie1<-boot.ci(x$boot.detail$x.new,matrix(x$boot.detail$ie1[,vari],nrow=length(x$boot.detail$x.new)),alpha,quantile)
    plot_ci(ie1,xlab=names(data$x)[data$dirx])
    
    if(!is.factor(data$x[,vari]))
   {if(full.model$distribution=="gaussian")
     suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else if(full.model$distribution=="coxph")
     suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    axis(1,at=data$x[,vari],labels=F)
    a<-marg.den(data$x[,data$dirx],data$x[,vari])
    #plot(a,xlab=names(data$x)[data$dirx],ylim=xlim,ylab="Mean.m")
    #lo <- loess(a[,1]~a[,2])
    #lines(lo$fitted[order(lo$x)], lo$x[order(lo$x)],lwd=1)
    scatter.smooth(a[,1],a[,2],family="gaussian",xlab=names(data$x)[data$dirx],ylim=xlim,ylab=paste("Mean",mname,sep="."))
    }
    
    #lo <- loess(data$x[,data$dirx]~data$x[,vari])
    #plot(data$x[,vari],data$x[,data$dirx],xlim=xlim,ylab=names(data$x)[data$dirx],xlab="")
    #lines(predict(lo), col='red', lwd=2)}
   else
   {if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
       suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
  }
}
else
{full.model=x$model$model
 data=x$data
 mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
 coef<-full.model$coefficients[names(full.model$coefficients)==vari] #plot the straight line instead of the loess line
   if(is.null(full.model$na.action))
    data1<-data$x[,vari]
  else
    data1<-data$x[-full.model$na.action,vari]
 if(data$binpred)
 {if(!is.factor(data$x[,vari]))
 {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  if(!x$model$Survival)
   b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values))
  else
    b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit)
  plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
  abline(a=mean(b[,2])-coef*mean(b[,1]),b=coef)
  #lo1 <- loess(b[,2]~b[,1])
  #lines(lo1$x[order(lo1$x)], lo1$fitted[order(lo1$x)], lwd=1)
    axis(1,at=data1,labels=F)
  overlapHist(a=data$x[,vari],b=data$x[,data$dirx],xlim=xlim,xname=names(data$x)[data$dirx])
  #temp<-data.frame(a=data$x[,vari],class= as.factor(data$x[,data$dirx]))
  #ggplot(temp, aes(a, fill = class)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
  #for (j in unique(data$x[,data$dirx]))
  #  hist(data$x[data$x[,data$dirx]==j,vari],ylab="Relative Freq",xlab="",xlim=xlim, freq=F,main=paste(names(data$x)[data$dirx],j,sep="="))
  #  j <- unique(data$x[,data$dirx])[1]
  #  p1<-hist(data$x[data$x[,data$dirx]==j,vari])
  #  j <- unique(data$x[,data$dirx])[2]
  #  p2<-hist(data$x[data$x[,data$dirx]==j,vari])
  #  plot( p1, col=rgb(0,0,1,1/4), ylab="Relative Freq",xlab="",xlim=xlim, freq=F,main=paste(names(data$x)[data$dirx],j,sep="="))  # first histogram
  #  plot( p2, col=rgb(1,0,0,1/4), xlim=xlim, add=T)  # second
 }
 else{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      if (!x$model$Survival)
        plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
      else
        plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
   #plot(data$x[,vari]~data$x[,data$dirx],ylab=mname,xlab=names(data$x)[data$dirx]) }
   temp1<-NULL
   for (j in unique(data$x[,data$dirx]))
     temp1<-c(temp1,prop.table(table(data$x[data$x[,data$dirx]==j,vari])))
   for (j in unique(data$x[,data$dirx]))
     barplot(prop.table(table(data$x[data$x[,data$dirx]==j,vari])),ylim=c(0,max(temp1)),
             ylab="Prop",sub=paste(names(data$x)[data$dirx],j,sep="="))}
 }
 else
 {par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
   ie1<-boot.ci(x$boot.detail$x.new,matrix(x$boot.detail$ie1[,vari],nrow=length(x$boot.detail$x.new)),alpha,quantile)
   plot_ci(ie1,xlab=names(data$x)[data$dirx])
   if(!is.factor(data$x[,vari]))
  {if(!x$model$Survival)
    b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values))
  else
    b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit)
#  if(is.null(full.model$na.action))
#    b<-marg.den(full.model$data[,vari],full.model$family$linkfun(full.model$fitted.values))
#   else
#    b<-marg.den(full.model$data[-full.model$na.action,vari],full.model$family$linkfun(full.model$fitted.values))
   plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
  # lo1 <- loess(b[,2]~b[,1])
  # lines(lo1$x[order(lo1$x)], lo1$fitted[order(lo1$x)], lwd=1)
   abline(a=mean(b[,2],na.rm=T)-coef*mean(b[,1],na.rm=T),b=coef)
   axis(1,at=data1,labels=F)
   a<-marg.den(data$x[,data$dirx],data$x[,vari])
   #plot(a,xlab=names(data$x)[data$dirx],ylim=xlim,ylab=paste("Mean",mname,sep="."))
   scatter.smooth(a[,1],a[,2],family="gaussian", xlab=names(data$x)[data$dirx],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
   #lo <- loess(a[,1]~a[,2])
   #lines(lo$fitted[order(lo$x)], lo$x[order(lo$x)],lwd=1)}
  else
  {if (!x$model$Survival)
     plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
   else  
     plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data$x[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
 }
}
par(op)
}
  



#plot on the med object
plot.med<-function(x,data,...,vari,xlim=range(data$x[,vari],na.rm=T))#data is the result from data.org
{marg.den<-function(x,y)
{y<-y[!is.na(x)]
x<-x[!is.na(x)]
z1<-unique(x)
z2<-rep(0,length(z1))
for (i in 1:length(z1))
  z2[i]<-mean(y[x==z1[i]],na.rm=T)
z3<-order(z1)
cbind(z1[z3],z2[z3])
}

overlapHist <- function(a, b,breaks=NULL, xlim=NULL, xname=NULL)
{a1<-a
b1<-b
a<-a[!is.na(a1) & !is.na(b1)]
b<-b[!is.na(a1) & !is.na(b1)]
j<-unique(b)
ahist<-hist(a[b==j[1]],plot=F)
dist = ahist$breaks[2]-ahist$breaks[1]
lb =min(ahist$breaks,na.rm = T)
ub=max(ahist$breaks,na.rm = T)
yl=max(ahist$density,na.rm = T)
for(i in j[-1])
{bhist<-hist(a[b==i],plot=F)
lb =min(lb,bhist$breaks,na.rm = T)
ub =max(ub,bhist$breaks,na.rm = T)
yl=max(yl,bhist$density,na.rm = T)
dist = min(dist,bhist$breaks[2]-bhist$breaks[1])
}
breaks=seq(lb,ub,dist)
if(is.null(xlim))
  xlim=c(lb,ub)
for (i in j)
  hist(a[b==i],ylab="Density",xlab="",breaks=breaks, 
       xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="="))
}


op <- par(no.readonly = TRUE) # the whole list of settable par's.
if (x$model[1]==T) 
{full.model=x$model$model
best.iter=x$model$best.iter
mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
if(data$binpred)
{if(!is.factor(data$x[,vari]))
{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
  else if(full.model$distribution=="coxph")
    suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
  else
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
  overlapHist(a=data$x[,vari],b=data$x[,data$dirx],xlim=xlim,xname=names(data$x)[data$dirx])
}
  else{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    temp1<-NULL
    for (j in unique(data$x[,data$dirx]))
      temp1<-c(temp1,prop.table(table(data$x[data$x[,data$dirx]==j,vari])))
    for (j in unique(data$x[,data$dirx]))
      barplot(prop.table(table(data$x[data$x[,data$dirx]==j,vari])),ylim=c(0,max(temp1,na.rm=T)),
              ylab="Prop",sub=paste(names(data$x)[data$dirx],j,sep="="))}
}
else
{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4)) #test
  temp2<-data$x[,data$dirx]
  temp3<-x$ie[,vari]
  temp.order=order(temp2)
 plot(temp2[temp.order], temp3[temp.order],type="l",
      xlab=names(data$x)[data$dirx],ylab=paste(c("IE of", vari),sep=""))
  
  if(!is.factor(data$x[,vari]))
  {if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    axis(1,at=data$x[,vari],labels=F)
    a<-marg.den(data$x[,data$dirx],data$x[,vari])
    scatter.smooth(a[,1],a[,2],family="gaussian",xlab=names(data$x)[data$dirx],ylim=xlim,ylab=paste("Mean",mname,sep="."))
  }
  
  else
  {if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
}
}
else
{full.model=x$model$model
mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
coef<-full.model$coefficients[names(full.model$coefficients)==vari] #plot the straight line instead of the loess line
if(is.null(full.model$na.action))
  data1<-data$x[,vari]
else
  data1<-data$x[-full.model$na.action,vari]
if(data$binpred)
{if(!is.factor(data$x[,vari]))
{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  if(!x$model$Survival)
    b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values))
  else
    b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit)
  plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
  abline(a=mean(b[,2])-coef*mean(b[,1]),b=coef)
  axis(1,at=data1,labels=F)
  overlapHist(a=data$x[,vari],b=data$x[,data$dirx],xlim=xlim,xname=names(data$x)[data$dirx])
}
  else{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if (!x$model$Survival)
      plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    else
      plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    temp1<-NULL
    for (j in unique(data$x[,data$dirx]))
      temp1<-c(temp1,prop.table(table(data$x[data$x[,data$dirx]==j,vari])))
    for (j in unique(data$x[,data$dirx]))
      barplot(prop.table(table(data$x[data$x[,data$dirx]==j,vari])),ylim=c(0,max(temp1)),
              ylab="Prop",sub=paste(names(data$x)[data$dirx],j,sep="="))}
}
else
{par(mfrow=c(3,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
  temp2<-data$x[,data$dirx]
  temp3<-x$ie[,vari]
  temp.order=order(temp2)
  plot(temp2[temp.order], temp3[temp.order],type="l",
       xlab=names(data$x)[data$dirx],ylab=paste(c("IE of", vari),sep=""))
  if(!is.factor(data$x[,vari]))
  {if(!x$model$Survival)
    b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values))
  else
    b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit)
  plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
  abline(a=mean(b[,2],na.rm=T)-coef*mean(b[,1],na.rm=T),b=coef)
  axis(1,at=data1,labels=F)
  a<-marg.den(data$x[,data$dirx],data$x[,vari])
  scatter.smooth(a[,1],a[,2],family="gaussian", xlab=names(data$x)[data$dirx],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
  else
  {if (!x$model$Survival)
    plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    else  
      plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data$x[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
    plot(data$x[,vari],data$x[,data$dirx],ylab=names(data$x)[data$dirx],xlab="")}
}
}
par(op)
}



################################################################################################
# the functions to make inferences on mediation effects using parallel calculation
med.par<-function(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm,contm = data$contm, 
                  catm = data$catm, jointm = data$jointm, allm = c(contm, catm), margin=1,
                  n=20,seed=sample(1:1000,1), nonlinear=F, df=1, nu=0.001,D=3,
                  distn=NULL,family1=NULL,refy=0,binpred=data$binpred,x.new=x,type=NULL,ncore=NULL)
{#for binary predictor
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, 
                     allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn="bernoulli",family1=binomial("logit"),
                     biny=F,refy=0,surv=F,type=NULL,ncore=NULL)
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {if(surv & !is.null(best.iter1))
      te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if (surv)
      te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    te
    }
    
    med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
    marg.m<-sample(marg.m)
    new1<-nom1
    new1[,med]<-marg.m[1:floor(n3/2)]
    new0<-nom0
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                 as.character(sample(nom0[,med],replace=T)))))
    else
      marg.m<-sample(c(sample(nom1[,med],replace=T),
                       sample(nom0[,med],replace=T)))
    marg.m<-sample(marg.m)}        #added for the new program
    else                                         #added for the new program
    {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
    new1<-nom1
    new0<-nom0
    if(length(med)==1)                                       #added for the new program, when there is only one mediator
    {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
    else                                                     #added for the new program
    {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
    marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
    dir.nom<-0
    for (i in levels(x[,med]))
    {new1<-nom1
    new1[1:dim(new1)[1],med]<-i
    new0<-nom0
    new0[1:dim(new0)[1],med]<-i
    p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
    if(surv & !is.null(best.iter1))
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
    else if(surv)
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
    else
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
    }
    dir.nom
    }
    
    #1.fit the model
    if (nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    registerDoParallel(cores=ncore)
    #2. repeat to get the mediation effect
    k<-1
    r<-foreach(k=1:n, .combine=rbind, .packages=c('gbm','survival','splines')) %dopar%
    {set.seed(seed+k)
      #2.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x1<-x[x[,dirx]==1,]
      x0<-x[x[,dirx]==0,]
      n3<-dim(x)[1]
      new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
      new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
      te<-te.binx(full.model,new1,new0,best.iter1,surv,type)
      denm<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
      #2.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm<-c(denm,med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type))
        }
      #2.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm<-c(denm,med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type))
        }
      #2.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm<-c(denm,med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type))
        }
      #2.5 get the indirect effects
      ie<-te-denm
      c(te,denm,ie)
    }
    stopImplicitCluster()
    
    #3. prepare for the store of results
    te<-r[,1]
    if(!is.null(jointm))
    {temp.col<-1+length(c(contm,catm))+jointm[[1]]
    temp.name<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
    else
    {temp.col<-1+length(c(contm,catm))
    temp.name<-c("de",names(x)[c(contm,catm)])}
    denm<-r[,2:(1+temp.col)]
    colnames(denm)<-temp.name
    ie<-r[,(2+temp.col):(ncol(r))]
    if(!is.null(jointm))
      dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
    else
      dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
    class(a)<-"med"
    return(a)
  }
  
  #for continous predictor
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                      biny=F,refy=0,x.new=x,surv=F,type=NULL,ncore=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    anymissing<-function(vec)
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row)
    {temp<-matrix(col,n.row)
    apply(temp,1,mean,na.rm=T)}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-x[,dirx]
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
    else
      models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
    {mult.norm<-function(mu,vari,n) 
    {  if (nrow(vari)!=ncol(vari)) 
    {result<-c("Error: Variance matrix is not square")}  
      else if (length(mu)!=nrow(vari)) 
      {result<-c("Error: length mu is not right!")}  
      else {   p<-length(mu)
      tmp1<-eigen(vari)$values
      tmp2<-eigen(vari)$vectors   
      result<-matrix(0,n,p)   
      for (i in 1:p)
      {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
      for (i in 1:n)
      {result[i,]<-tmp2%*%result[i,]+mu}
      }  
      result
    }
    
    match.margin<-function(vec)   
    {range1<-vec[1:2]
    vec1<-vec[-(1:2)]
    range2<-range(vec1,na.rm=T)
    vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
    vec1
    }
    
    gen.mult<-function(vec)
    {if(sum(is.na(vec))>0)
      return(rep(NA,length(vec)))
      else{ 
        l<-1-sum(vec)
        l<-ifelse(l<0,0,l)
        return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
    }
    
    means<-NULL
    z<-x1[,dirx]
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
    if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
    {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
    sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
    else{
      sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
      
      range.means<-apply(means,2,range,na.rm=T)
      
      sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
    }
    sim.m2<-data.frame(sim.m2)
    n<-dim(sim.m2)[1]
    if(!is.null(binm))
      for (i in 1:length(binm))
        sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
    
    if(!is.null(catm))
    {j<-length(binm)+1
    for (i in 2:(catm$n+1))
    {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
    sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
    j<-j+length(catm[[i]])}
    }
    
    x1[,c(binm1,contm)]<-sim.m2
    
    x1
    }
    
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
    x<-x[nonmissing,]
    y<-y[nonmissing]
    nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
    x.new<-x.new[nonmissing1,]
    
    #1.fit the model
    if(nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
    }
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    #2. prepare for the store of results
    n.new<-nrow(x.new)
    te<-rep(0,n.new)
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
    #te1<-NULL
    #x1<-x.new
    #x1[,dirx]<-x[,dirx]+margin
    #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
    
    n1<-dim(x)[1]
    #denm1<-NULL
    
    #4. repeat to get the mediation effect
    registerDoParallel(cores=ncore)
   k<-1
    r<-foreach(k=1:n, .combine=rbind, .packages=c('gbm','survival','splines')) %dopar%
    {set.seed(seed+k)
      new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
      x1<-new0
      x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
      if(surv & !is.null(best.iter1))
      {ybar0<-predict(full.model,new0,best.iter1,type=type)
      denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
      else if(surv)
      {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
      denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
      else
      {ybar0<-predict(full.model,new0,best.iter1)
      denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
      new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
      #4.1 get the te
      if(surv & !is.null(best.iter1))
        te1<-(predict(full.model,new1,best.iter1,type=type)-ybar0)/margin
      else if(surv)
        te1<-(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin
      else
        te1<-(predict(full.model,new1,best.iter1)-ybar0)/margin
      
      #4.2 mediation effect from the single mediator
      if (!is.null(listm$single))
        for (i in 1:length(listm$single))
        {new1.nm<-new1
        new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
        if(surv & !is.null(best.iter1))
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
        else if(surv)
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
        }
      
      #4.3.mediation effect from the joint mediator
      if (!is.null(listm$multi))
        for (i in 2:(listm$multi[[1]]+1))
        {new1.nm<-new1
        new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
        if(surv & !is.null(best.iter1))
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
        else if(surv)
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
        }
      cbind(te1,denm2)
    } 
    stopImplicitCluster()
    
    #4.4 get the indirect effects
    te1<-matrix(r[,1],n.new)
    denm1<-r[,-1]
    denm<-apply(denm1,2,col_mean,n.new)
    te<-apply(te1,1,mean)
    
    ie<-te-denm
    if(!is.null(listm$multi))
      colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(denm)<-c("de",names(x)[listm$single])
    if(!is.null(listm$multi))
      colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(ie)<-c("all",names(x)[listm$single])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new)
    class(a)<-"med"
    return(a)
  }
  
  surv=F
  if(class(y)=="Surv")
  {surv=T
  biny=F
  if(is.null(distn))
    distn = "coxph"
  if(is.null(type) & nonlinear)
    type="response"
  else if (is.null(type))
    type="risk"}
  else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
  {biny=T
  if(is.null(family1))
    family1 = binomial("logit")
  if(is.null(distn))
    distn = "bernoulli"}
  else
  {biny=F
  if(is.null(family1))
    family1 = gaussian(link = "identity")
  if(is.null(distn))
    distn = "gaussian"}
  
  if(binpred)
    a<-med.binx(data=data, x=x, y=y, dirx=dirx, contm = contm, 
                catm=catm, jointm=jointm, allm=allm, n=n,seed=seed,
                nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,
                biny=biny,refy=refy,surv=surv,type=type)
  else
    a<-med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                 catm=catm, jointm=jointm, margin=margin, n=n, seed=seed, 
                 nonlinear=nonlinear, df=df, nu=nu,D=D, distn=distn, 
                 family1=family1,biny=biny,refy=refy,x.new=x.new,surv=surv,type=type)
  return(a)
}


boot.med.par<-function(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,catm=data$catm,
                       jointm=data$jointm,margin=1,n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,
                       D=3,distn=NULL,family1=NULL,n2=50,weight=rep(1,nrow(x)),refy=NULL,
                       x.new=x,binpred=data$binpred,type=NULL,ncore=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         weight=rep(1,nrow(x)),biny=F,refy=0,surv,type,ncore=NULL)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, 
                     allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn="bernoulli",family1=binomial("logit"),
                     biny=F,refy=0,surv=F,type=NULL)
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {if(surv & !is.null(best.iter1))
      te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if (surv)
      te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    te
    }
    
    med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
    marg.m<-sample(marg.m)
    new1<-nom1
    new1[,med]<-marg.m[1:floor(n3/2)]
    new0<-nom0
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                 as.character(sample(nom0[,med],replace=T)))))
    else
      marg.m<-sample(c(sample(nom1[,med],replace=T),
                       sample(nom0[,med],replace=T)))
    marg.m<-sample(marg.m)}        #added for the new program
    else                                         #added for the new program
    {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
    new1<-nom1
    new0<-nom0
    if(length(med)==1)                                       #added for the new program, when there is only one mediator
    {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
    else                                                     #added for the new program
    {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
    marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
    dir.nom<-0
    for (i in levels(x[,med]))
    {new1<-nom1
    new1[1:dim(new1)[1],med]<-i
    new0<-nom0
    new0[1:dim(new0)[1],med]<-i
    p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
    if(surv & !is.null(best.iter1))
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
    else if(surv)
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
    else
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
    }
    dir.nom
    }
    
    #1.fit the model
    if (nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    #2. prepare for the store of results
    set.seed(seed)
    te<-rep(0,n)
    if(!is.null(jointm))
    {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
    else
    {denm<-matrix(0,n,1+length(c(contm,catm)))
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x1<-x[x[,dirx]==1,]
      x0<-x[x[,dirx]==0,]
      n3<-dim(x)[1]
      new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
      new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
      te[k]<-te.binx(full.model,new1,new0,best.iter1,surv,type)
      denm[k,1]<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type)
        j<-j+1}
      
      #3.5 get the indirect effects
      ie[k,]<-te[k]-denm[k,]
      if(!is.null(jointm))
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
    }
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  if(is.null(jointm))
  {ie<-matrix(0,n2+1,1+length(c(contm,catm)))
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])}
  else
  {ie<-matrix(0,n2+1,1+length(c(contm,catm))+jointm[[1]])
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  temp<-med.binx(data=F,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1]<-mean(temp$te)
  de[1]<-mean(temp$denm[,1])
  ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  model<-temp$model
  best.iter<-temp$best.iter
  registerDoParallel(cores=ncore)
  r<-foreach(i=1:n2,.combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T,prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.binx(data=F,x1,y1,dirx,contm,catm,jointm,allm,n,seed+i,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te1<-mean(temp$te)
  de1<-mean(temp$denm[,1])
  ie1<-apply(temp$ie,2,mean)  #first row is the estimated value
  c(te1,de1,ie1)  
  }
  
  te[2:(n2+1)]<-r[,1]
  de[2:(n2+1)]<-r[,2]
  ie[2:(n2+1),]<-r[,-c(1,2)]
  
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         weight=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type,ncore=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                      biny=F,refy=0,x.new=x,surv=F,type=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    anymissing<-function(vec)
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row)
    {temp<-matrix(col,n.row)
    apply(temp,1,mean,na.rm=T)}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-x[,dirx]
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
    else
      models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
    {mult.norm<-function(mu,vari,n) 
    {  if (nrow(vari)!=ncol(vari)) 
    {result<-c("Error: Variance matrix is not square")}  
      else if (length(mu)!=nrow(vari)) 
      {result<-c("Error: length mu is not right!")}  
      else {   p<-length(mu)
      tmp1<-eigen(vari)$values
      tmp2<-eigen(vari)$vectors   
      result<-matrix(0,n,p)   
      for (i in 1:p)
      {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
      for (i in 1:n)
      {result[i,]<-tmp2%*%result[i,]+mu}
      }  
      result
    }
    
    match.margin<-function(vec)   
    {range1<-vec[1:2]
    vec1<-vec[-(1:2)]
    range2<-range(vec1,na.rm=T)
    vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
    vec1
    }
    
    gen.mult<-function(vec)
    {if(sum(is.na(vec))>0)
      return(rep(NA,length(vec)))
      else{ 
        l<-1-sum(vec)
        l<-ifelse(l<0,0,l)
        return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
    }
    
    means<-NULL
    z<-x1[,dirx]
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
    if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
    {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
    sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
    else{
      sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
      
      range.means<-apply(means,2,range,na.rm=T)
      
      sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
    }
    sim.m2<-data.frame(sim.m2)
    n<-dim(sim.m2)[1]
    if(!is.null(binm))
      for (i in 1:length(binm))
        sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
    
    if(!is.null(catm))
    {j<-length(binm)+1
    for (i in 2:(catm$n+1))
    {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
    sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
    j<-j+length(catm[[i]])}
    }
    
    x1[,c(binm1,contm)]<-sim.m2
    
    x1
    }
    
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
    x<-x[nonmissing,]
    y<-y[nonmissing]
    nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
    x.new<-x.new[nonmissing1,]
    
    #1.fit the model
    if(nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
    }
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    te<-rep(0,n.new)
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
    te1<-NULL
    #x1<-x.new
    #x1[,dirx]<-x[,dirx]+margin
    #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
    
    n1<-dim(x)[1]
    denm1<-NULL
    
    #4. repeat to get the mediation effect
    for (k in 1:n)
    {new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
    x1<-new0
    x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
    if(surv & !is.null(best.iter1))
    {ybar0<-predict(full.model,new0,best.iter1,type=type)
    denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
    else if(surv)
    {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
    denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
    else
    {ybar0<-predict(full.model,new0,best.iter1)
    denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
    new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
    #4.1 get the te
    if(surv & !is.null(best.iter1))
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      te1<-cbind(te1, (predict(full.model,new1,best.iter1)-ybar0)/margin)
    
    #4.2 mediation effect from the single mediator
    if (!is.null(listm$single))
      for (i in 1:length(listm$single))
      {new1.nm<-new1
      new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
      for (i in 2:(listm$multi[[1]]+1))
      {new1.nm<-new1
      new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    denm1<-rbind(denm1,denm2)
    } 
    
    #4.4 get the indirect effects
    denm<-apply(denm1,2,col_mean,n.new)
    te<-apply(te1,1,mean)
    ie<-te-denm
    if(!is.null(listm$multi))
      colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(denm)<-c("de",names(x)[listm$single])
    if(!is.null(listm$multi))
      colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(ie)<-c("all",names(x)[listm$single])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new)
    class(a)<-"med"
    return(a)
  }
  
  col_mean1<-function(col,n.row)
  {temp<-matrix(col,n.row)
  apply(temp,2,mean,na.rm=T)}
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(binm))
    binm<-unlist(sapply(binm,grep,xnames))
  if(!is.null(catm))
    for (i in 2:length(catm))
      if(is.character(catm[[i]]))
        catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  if(is.null(catm))
  {multi=jointm
  name1<-NULL                       #added in the new program
  if (!is.null(multi))              #added in the new program, in case that multi is NULL
    name1<-paste("j",1:multi[[1]],sep="")}
  else if(is.null(jointm))
  {multi=catm
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])}
  else {temp1<-catm
  temp2<-jointm
  temp1[[1]]=catm[[1]]+jointm[[1]]
  temp2[[1]]<-NULL
  multi=append(temp1,temp2)
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])
  name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
  listm=list(single=c(contm,binm),multi=multi)
  
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2+1,1+length(listm$single)+mul)   #added in the new program
  if(!is.null(listm$multi))
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single],name1)
  else 
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single])
  
  temp<-med.contx(data=F,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,surv,type)
  x.new<-temp$x.new
  te[1]<-mean(temp$te,na.rm=T)
  de[1]<-mean(temp$denm[,1],na.rm=T) 
  ie[1,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  model<-temp$model
  registerDoParallel(cores=ncore)
  l<-1
  r<-foreach(l=1:n2, .combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T, prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.contx(data=F,x1,y1,dirx,binm,contm,catm,jointm, margin,n,seed+l,nonlinear,df,nu,D,distn,family1,biny,refy,x.new,surv,type) #added to the new codel, change the seed to make different results
  #te[1+l]<-mean(temp$te,na.rm=T)
  #de[1+l]<-mean(temp$denm[,1],na.rm=T)
  #ie[1+l,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  te1<-temp$te
  de1<-temp$denm[,1]
  ie1<-temp$ie
  cbind(te1,de1,ie1)
  }
  te1<-matrix(r[,1],nrow(x.new))
  de1<-matrix(r[,2],nrow(x.new))
  ie1<-r[,-c(1,2)]
  te[2:(n2+1)]<-apply(te1,2,mean)
  de[2:(n2+1)]<-apply(de1,2,mean)
  ie[2:(n2+1),]<-apply(ie1,2,col_mean1,nrow(x.new))
  colnames(ie1)<-colnames(ie)
  
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(x.new=x.new[,dirx],te1=te1,de1=de1,ie1=ie1))
  class(a)<-"mma"
  return(a)
}

surv=F
if(class(y)=="Surv")
{surv=T
biny=F
if(is.null(distn))
  distn = "coxph"
if(is.null(type) & nonlinear)
  type="response"
else if (is.null(type))
  type="risk"}
else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
{biny=T
if(is.null(family1))
  family1 = binomial("logit")
if(is.null(distn))
  distn = "bernoulli"
if(!is.null(refy))
  y<-ifelse(y==refy,0,1)
else
  y<-ifelse(as.factor(y)==levels(as.factor(y))[1],0,1)
}
else
{biny=F
if(is.null(family1))
  family1 = gaussian(link = "identity")
if(is.null(distn))
  distn = "gaussian"
}

if(binpred)
  a<-boot.med.binx(data=data,x=x, y=y,dirx=dirx,contm=contm,catm=catm,
                   jointm=jointm,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,
                   D=D,distn=distn,family1=family1,
                   weight=weight,biny=biny,refy=0,surv,type,ncore=ncore)
else
  a<-boot.med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                    catm=catm, jointm=jointm, margin = margin, n = n, seed = seed, 
                    nonlinear = nonlinear, df = df, nu = nu, D = D, distn = distn, 
                    family1 = family1, n2 = n2,weight=weight,
                    biny=biny,refy=0,x.new=x.new,surv,type,ncore=ncore)

return(a)
}


mma.par<-function(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,
                  catmed=NULL,catref=NULL,jointm=NULL,refy=NULL,
                  predref=NULL,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
                  nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=NULL,n2=50,
                  weight=rep(1,nrow(x)),x.new=NULL,type=NULL,ncore=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         weight=rep(1,nrow(x)),biny=F,refy=0,surv,type,ncore=NULL)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, 
                     allm = c(contm, catm), n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn="bernoulli",family1=binomial("logit"),
                     biny=F,refy=0,surv=F,type=NULL)
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {if(surv & !is.null(best.iter1))
      te<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if (surv)
      te<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      te<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    te
    }
    
    med.binx.contm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m<-c(nom1[sample(1:dim(nom1)[1],replace=T),med],nom0[sample(1:dim(nom0)[1],replace=T),med])
    marg.m<-sample(marg.m)
    new1<-nom1
    new1[,med]<-marg.m[1:floor(n3/2)]
    new0<-nom0
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-sample(as.factor(c(as.character(sample(nom1[,med],replace=T)),
                                 as.character(sample(nom0[,med],replace=T)))))
    else
      marg.m<-sample(c(sample(nom1[,med],replace=T),
                       sample(nom0[,med],replace=T)))
    marg.m<-sample(marg.m)}        #added for the new program
    else                                         #added for the new program
    {marg.m<-rbind(nom1[sample(1:nrow(nom1),replace=T),med],nom0[sample(1:nrow(x0),replace=T),med])
    marg.m<-marg.m[sample(2*floor(n3/2)),]  }     
    new1<-nom1
    new0<-nom0
    if(length(med)==1)                                       #added for the new program, when there is only one mediator
    {new1[,med]<-marg.m[1:floor(n3/2)]                     #added for the new program 
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2))]}  #added for the new program
    else                                                     #added for the new program
    {new1[,med]<-marg.m[1:floor(n3/2),]
    new0[,med]<-marg.m[(floor(n3/2)+1):(2*floor(n3/2)),]}
    if(surv & !is.null(best.iter1))
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T)
    else if(surv)
      dir.nom<-mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom<-mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T)
    dir.nom
    }
    
    med.binx.catm<-function(full.model,x,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-dim(x)[1]
    marg.m1<-nom1[sample(dim(nom1)[1],floor(n3/2),replace=T),med]
    marg.m2<-nom0[sample(dim(nom0)[1],floor(n3/2),replace=T),med]
    dir.nom<-0
    for (i in levels(x[,med]))
    {new1<-nom1
    new1[1:dim(new1)[1],med]<-i
    new0<-nom0
    new0[1:dim(new0)[1],med]<-i
    p<-0.5*mean(marg.m1==i,na.rm=T)+ 0.5*mean(marg.m2==i,na.rm=T)
    if(surv & !is.null(best.iter1))
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type),na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type),na.rm=T))
    else if(surv)
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit,na.rm=T))
    else
      dir.nom<-dir.nom+p*(mean(predict(full.model,new1,best.iter1),na.rm=T)- mean(predict(full.model,new0,best.iter1),na.rm=T))
    }
    dir.nom
    }
    
    #1.fit the model
    if (nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}}
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    #2. prepare for the store of results
    set.seed(seed)
    te<-rep(0,n)
    if(!is.null(jointm))
    {denm<-matrix(0,n,1+length(c(contm,catm))+jointm[[1]])
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
    else
    {denm<-matrix(0,n,1+length(c(contm,catm)))
    dimnames(denm)[[2]]<-c("de",names(x)[c(contm,catm)])}
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x1<-x[x[,dirx]==1,]
      x0<-x[x[,dirx]==0,]
      n3<-dim(x)[1]
      new1<-x1[sample(1:dim(x1)[1],floor(n3/2),replace=T),]
      new0<-x0[sample(1:dim(x0)[1],floor(n3/2),replace=T),]
      te[k]<-te.binx(full.model,new1,new0,best.iter1,surv,type)
      denm[k,1]<-med.binx.jointm(full.model,x,new1,new0,allm,best.iter1,surv,type)
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.contm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.catm(full.model,x,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[k,j]<-med.binx.jointm(full.model,x,new1,new0,jointm[[i+1]],best.iter1,surv,type)
        j<-j+1}
      
      #3.5 get the indirect effects
      ie[k,]<-te[k]-denm[k,]
      if(!is.null(jointm))
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])
    }
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1))
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  if(is.null(jointm))
  {ie<-matrix(0,n2+1,1+length(c(contm,catm)))
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)])}
  else
  {ie<-matrix(0,n2+1,1+length(c(contm,catm))+jointm[[1]])
  dimnames(ie)[[2]]<-c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))}
  temp<-med.binx(data=F,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te[1]<-mean(temp$te)
  de[1]<-mean(temp$denm[,1])
  ie[1,]<-apply(temp$ie,2,mean)  #first row is the estimated value
  model<-temp$model
  best.iter<-temp$best.iter
  registerDoParallel(cores=ncore)
  r<-foreach(i=1:n2,.combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T,prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.binx(data=F,x1,y1,dirx,contm,catm,jointm,allm,n,seed+i,nonlinear,nu,D,distn,family1,biny,refy,surv,type)
  te1<-mean(temp$te)
  de1<-mean(temp$denm[,1])
  ie1<-apply(temp$ie,2,mean)  #first row is the estimated value
  c(te1,de1,ie1)  
  }
  
  te[2:(n2+1)]<-r[,1]
  de[2:(n2+1)]<-r[,2]
  ie[2:(n2+1),]<-r[,-c(1,2)]
  
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         weight=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type,ncore=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",family1=gaussian(link="identity"),
                      biny=F,refy=0,x.new=x,surv=F,type=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    if(is.character(dirx))
      pred<-grep(dirx,xnames)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    if(biny)                     #recode y if y is binary
      y<-ifelse(y==refy,0,1)
    
    x<-x[!is.na(y),]             #delete nas in y for mart
    y<-y[!is.na(y)]
    
    anymissing<-function(vec)
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row)
    {temp<-matrix(col,n.row)
    apply(temp,1,mean,na.rm=T)}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df)
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-x[,dirx]
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=binomial(link = "logit"))
    else
      models[[j]]<-glm(x[,i]~z,family=binomial(link = "logit"))
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    {if(nonlinear)
      models[[j]]<-glm(x[,i]~ns(z,df=df),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(x[,i]~z,family=gaussian(link="identity"))
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm)
    {mult.norm<-function(mu,vari,n) 
    {  if (nrow(vari)!=ncol(vari)) 
    {result<-c("Error: Variance matrix is not square")}  
      else if (length(mu)!=nrow(vari)) 
      {result<-c("Error: length mu is not right!")}  
      else {   p<-length(mu)
      tmp1<-eigen(vari)$values
      tmp2<-eigen(vari)$vectors   
      result<-matrix(0,n,p)   
      for (i in 1:p)
      {result[,i]<-rnorm(n,mean=0,sd=sqrt(tmp1[i]))}   
      for (i in 1:n)
      {result[i,]<-tmp2%*%result[i,]+mu}
      }  
      result
    }
    
    match.margin<-function(vec)   
    {range1<-vec[1:2]
    vec1<-vec[-(1:2)]
    range2<-range(vec1,na.rm=T)
    vec1<-range1[1]+diff(range1)/diff(range2)*(vec1-range2[1])
    vec1
    }
    
    gen.mult<-function(vec)
    {if(sum(is.na(vec))>0)
      return(rep(NA,length(vec)))
      else{ 
        l<-1-sum(vec)
        l<-ifelse(l<0,0,l)
        return(rmultinom(1,size=1,prob=c(l,vec))[-1])}
    }
    
    means<-NULL
    z<-x1[,dirx]
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z=z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z=z)))
    if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
    {sim.m<-rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat))     #added in the new program
    sim.m2<-match.margin(c(range(means,na.rm=T),sim.m))}                          #added in the new program   
    else{
      sim.m<-t(apply(means,1,mult.norm,vari=distmgivenx$varmat,n=1))
      
      range.means<-apply(means,2,range,na.rm=T)
      
      sim.m2<-apply(rbind(range.means,sim.m),2,match.margin)    #to make the simulate fit the means' ranges
    }
    sim.m2<-data.frame(sim.m2)
    n<-dim(sim.m2)[1]
    if(!is.null(binm))
      for (i in 1:length(binm))
        sim.m2[,i]<-rbinom(n,size=1,prob=sim.m2[,i])
    
    if(!is.null(catm))
    {j<-length(binm)+1
    for (i in 2:(catm$n+1))
    {a<-sim.m2[,j:(j+length(catm[[i]])-1)]
    sim.m2[,j:(j+length(catm[[i]])-1)]<-t(apply(a,1,gen.mult))
    j<-j+length(catm[[i]])}
    }
    
    x1[,c(binm1,contm)]<-sim.m2
    
    x1
    }
    
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    nonmissing<-apply(cbind(y,x[,c(dirx,listm$single)]),1,anymissing)
    x<-x[nonmissing,]
    y<-y[nonmissing]
    nonmissing1<-apply(x.new[,c(dirx,listm$single)],1,anymissing)
    x.new<-x.new[nonmissing1,]
    
    #1.fit the model
    if(nonlinear)
    {full.model<-suppressWarnings(gbm.fit(x,y, n.trees=200, interaction.depth=D, shrinkage=nu,
                                          distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
    best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))         
    while(full.model$n.trees-best.iter1<30){
      full.model<-gbm.more(full.model, 50)           # do another 50 iterations
      best.iter1<-suppressWarnings(gbm.perf(full.model,plot.it=FALSE,method="OOB"))}
    }
    else
    {if(surv)
      full.model<-coxph(y~., data=x)
    else
      full.model<-glm(y~., data=x, family=family1)
    best.iter1=NULL}
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    te<-rep(0,n.new)
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,dirx,binm,contm,catm,nonlinear,df)
    te1<-NULL
    #x1<-x.new
    #x1[,dirx]<-x[,dirx]+margin
    #ybar0<-mean(predict(full.model,x,best.iter1),na.rm=T)
    
    n1<-dim(x)[1]
    denm1<-NULL
    
    #4. repeat to get the mediation effect
    for (k in 1:n)
    {new0<-sim.xm(distmgivenx,x.new,dirx,binm,contm,catm) #draw ms conditional on x.new
    x1<-new0
    x1[,dirx]<-new0[,dirx]+margin #assume can change x without changing others
    if(surv & !is.null(best.iter1))
    {ybar0<-predict(full.model,new0,best.iter1,type=type)
    denm2<-(predict(full.model,x1,best.iter1,type=type)-ybar0)/margin}
    else if(surv)
    {ybar0<-predict(full.model,new0,best.iter1,type=type,se.fit=TRUE)$fit
    denm2<-(predict(full.model,x1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin}
    else
    {ybar0<-predict(full.model,new0,best.iter1)
    denm2<-(predict(full.model,x1,best.iter1)-ybar0)/margin}
    new1<-sim.xm(distmgivenx,x1,dirx,binm,contm,catm)  #draw from the conditional distribution of m given x
    #4.1 get the te
    if(surv & !is.null(best.iter1))
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type)-ybar0)/margin)
    else if(surv)
      te1<-cbind(te1, (predict(full.model,new1,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
    else
      te1<-cbind(te1, (predict(full.model,new1,best.iter1)-ybar0)/margin)
    
    #4.2 mediation effect from the single mediator
    if (!is.null(listm$single))
      for (i in 1:length(listm$single))
      {new1.nm<-new1
      new1.nm[,listm$single[i]]<-new0[,listm$single[i]]    #draw m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
      for (i in 2:(listm$multi[[1]]+1))
      {new1.nm<-new1
      new1.nm[,listm$multi[[i]]]<-new0[,listm$multi[[i]]]    #draw joint m from its original distribution
      if(surv & !is.null(best.iter1))
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type)-ybar0)/margin)
      else if(surv)
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1,type=type,se.fit=TRUE)$fit-ybar0)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model,new1.nm,best.iter1)-ybar0)/margin)
      }
    denm1<-rbind(denm1,denm2)
    } 
    
    #4.4 get the indirect effects
    denm<-apply(denm1,2,col_mean,n.new)
    te<-apply(te1,1,mean)
    ie<-te-denm
    if(!is.null(listm$multi))
      colnames(denm)<-c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(denm)<-c("de",names(x)[listm$single])
    if(!is.null(listm$multi))
      colnames(ie)<-c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep=""))
    else 
      colnames(ie)<-c("all",names(x)[listm$single])
    
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new)
    class(a)<-"med"
    return(a)
  }
  
  col_mean1<-function(col,n.row)
  {temp<-matrix(col,n.row)
   apply(temp,2,mean,na.rm=T)}
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  if(is.character(dirx))
    pred<-grep(dirx,xnames)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(binm))
    binm<-unlist(sapply(binm,grep,xnames))
  if(!is.null(catm))
    for (i in 2:length(catm))
      if(is.character(catm[[i]]))
        catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  if(is.null(catm))
  {multi=jointm
  name1<-NULL                       #added in the new program
  if (!is.null(multi))              #added in the new program, in case that multi is NULL
    name1<-paste("j",1:multi[[1]],sep="")}
  else if(is.null(jointm))
  {multi=catm
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])}
  else {temp1<-catm
  temp2<-jointm
  temp1[[1]]=catm[[1]]+jointm[[1]]
  temp2[[1]]<-NULL
  multi=append(temp1,temp2)
  name1<-NULL
  for (i in 2:(catm[[1]]+1))
    name1<-c(name1,names(x)[multi[[i]][1]])
  name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
  listm=list(single=c(contm,binm),multi=multi)
  
  te<-rep(0,n2+1)
  de<-rep(0,n2+1)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2+1,1+length(listm$single)+mul)   #added in the new program
  if(!is.null(listm$multi))
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single],name1)
  else 
    dimnames(ie)[[2]]<-c("all",names(x)[listm$single])
  
  temp<-med.contx(data=F,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,surv,type)
  x.new=temp$x.new
  te[1]<-mean(temp$te,na.rm=T)
  de[1]<-mean(temp$denm[,1],na.rm=T) 
  ie[1,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  te1<-NULL                      #to store the mediation effects on predictor
  de1<-NULL
  ie1<-NULL
  model<-temp$model
  registerDoParallel(cores=ncore)
  l<-1
  r<-foreach(l=1:n2, .combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T, prob=weight)
  x1<-x[boots,]
  y1<-y[boots]
  temp<-med.contx(data=F,x1,y1,dirx,binm,contm,catm,jointm, margin,n,seed+l,nonlinear,df,nu,D,distn,family1,biny,refy,x.new,surv,type) #added to the new codel, change the seed to make different results
  #te[1+l]<-mean(temp$te,na.rm=T)
  #de[1+l]<-mean(temp$denm[,1],na.rm=T)
  #ie[1+l,]<-apply(temp$ie,2,mean,na.rm=T)  #first row is the estimated value
  te1<-temp$te
  de1<-temp$denm[,1]
  ie1<-temp$ie
  cbind(te1,de1,ie1)
  }
  te1<-matrix(r[,1],nrow(x.new))
  de1<-matrix(r[,2],nrow(x.new))
  ie1<-r[,-c(1,2)]
  te[2:(n2+1)]<-apply(te1,2,mean)
  de[2:(n2+1)]<-apply(de1,2,mean)
  ie[2:(n2+1),]<-apply(ie1,2,col_mean1,nrow(x.new))
  colnames(ie1)<-colnames(ie)
  
  a<-list(estimation=list(ie=ie[1,],te=te[1],de=de[1]),bootsresults=list(ie=ie[-1,],te=te[-1],de=de[-1]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(x.new=x.new[,dirx],te1=te1,de1=de1,ie1=ie1))
  class(a)<-"mma"
  return(a)
}


surv=F
if(class(y)=="Surv")
{surv=T
biny=F
if(is.null(distn))
  distn = "coxph"
if(is.null(type) & nonlinear)
  type="response"
else if (is.null(type))
  type="risk"}
else if(is.character(y) | is.factor(y) | nlevels(as.factor(y))==2)
{biny=T
if(is.null(family1))
  family1 = binomial("logit")
if(is.null(distn))
  distn = "bernoulli"
if(!is.null(refy))
  y<-ifelse(y==refy,0,1)
else
  y<-ifelse(as.factor(y)==levels(as.factor(y))[1],0,1)
}
else
{biny=F
if(is.null(family1))
  family1 = gaussian(link = "identity")
if(is.null(distn))
  distn = "gaussian"
}

data<-data.org(x=x,y=y,pred=pred,mediator=mediator,contmed=contmed,binmed=binmed,
               binref=binref,catmed=catmed,catref=catref,jointm=jointm,refy=refy,family1=family1,
               predref=predref,alpha=alpha,alpha2=alpha2)
binpred<-data$binpred

if(binpred) 
{result<-boot.med.binx(data=data,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,weight=weight,biny=biny,refy=0,surv=surv,type=type,ncore=ncore)
}
else
{if(is.null(x.new))
  result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                         D=D,distn=distn,family1=family1,n2=n2,weight=weight,biny=biny,refy=0,surv=surv,type=type,ncore=ncore)
else
  result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                         D=D,distn=distn,family1=family1,n2=n2,weight=weight,biny=biny,refy=0, x.new=x.new,surv=surv,type=type,ncore=ncore)

}
result
}

  
