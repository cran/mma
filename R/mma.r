#to organize data
data.org<-function(x,y,pred,mediator=NULL,contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
                   catref=NULL,jointm=NULL,refy=rep(NA,ncol(data.frame(y))), 
                   family1=as.list(rep(NA,ncol(data.frame(y)))),
                   predref=NULL,alpha=0.1,alpha2=0.1,testtype=1, w=NULL)
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
 {a<-factor(x[,i])
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

#y2 is the response dataframe
y2<-data.frame(y)                     #consider multivarite or multicategorical responses
ny<-ncol(y2)
y_type<-rep(4,ny)                     #1 is continuous, 2 is binary, 3 is multi-categorical, 4 is survival
for (i in 1:ny)
{if(class(y2[,i])!="Surv")
  if(nlevels(droplevels(as.factor(y2[,i])))==2)   
  {y_type[i]<-2
   if(is.na(family1[[i]]))
    family1[[i]] = binomial("logit")
   if(!is.na(refy[i]))
    y2[,i]<-ifelse(y2[,i]==refy[i],0,1)
  else
   {refy[i]<-levels(droplevels(as.factor(y2[,i])))[1]
    y2[,i]<-ifelse(as.factor(y2[,i])==refy[i],0,1)}
 }
 else if(is.character(y2[,i]) | is.factor(y2[,i]))
 {y_type[i]<-3
  y2[,i]<-droplevels(y2[,i])  #drop those empty levles
  if(is.na(refy[i]))
    refy[i]<-levels(as.factor(y2[,i]))[1]
 }
 else
 {y_type[i]=1
  if(is.na(family1[[i]]))
   family1[[i]] = gaussian(link = "identity")
 }
}

if(sum(y_type==3)>0) #transfer the multicategorical type response
{temp1<-(1:length(y_type))[y_type==3]
 temp2<-cattobin(y2,temp1,refy[temp1])
 y2<-data.frame(temp2$x)
 y_type<-y_type[-temp1]
 y_type<-c(y_type,rep(2,ncol(y2)-length(y_type)))
 family1[[temp1]]<-NULL
 family1<-append(family1,rep(list(binomial("logit")),ncol(y2)-length(family1)))
}

xnames<-colnames(x)

#predictors have to be one categorical or all continuous, pred is the exposure vector/matrix
pred1<-data.frame(pred)
binpred=T
for (i in 1:ncol(pred1))
 if(nlevels(as.factor(pred1[,i]))==2)
 {if(!is.null(predref))
      pred1[,i]<-as.factor(ifelse(pred1[,i]==predref,0,1))
  else
      {pred1[,i]<-as.factor(pred1[,i])
       pred1[,i]<-as.factor(ifelse(pred1[,i]==levels(pred1[,i])[1],0,1))}
 }
 else if(is.character(pred1[,i]) | is.factor(pred1[,i]))
 {if(!is.null(predref))
   pred1<-cattobin(data.frame(pred1[,i]),1,predref)$x
  else
   pred1<-cattobin(data.frame(pred1[,i]),1,levels(as.factor(pred1[,i]))[1])$x
 }
else
  binpred=F
pred<-data.frame(pred1)

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

if(!is.null(binmed))  #revise: binref does not have to be null or full, & is.null(binref)
 {j<-1
  for(i in binmed)
  {x[,i]=as.factor(x[,i])
   if(is.null(binref))
   binref[i]=levels(x[,i])[1]
   else if(is.na(binref[i]))
     binref[i]=levels(x[,i])[1]
   j<-j+1}}

if(!is.null(catmed))   #revise: catref does not have to be null or full, & is.null(catref)
 {j<-1  
  for(i in catmed)
  {x[,i]=as.factor(x[,i])
    if(is.null(catref))
      catref[i]=levels(x[,i])[1]
    else if(is.na(catref[j]))
      catref[i]=levels(x[,i])[1]
   j<-j+1}}

if(!is.null(mediator))   #for all potential mediators
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

if(!is.null(jointm))  #mediators that should be jointly considered are forced in as mediators
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

#newx1<-cbind(newx1,pred) #add the predictor(s) back
#delete variables that are not significant
fullmodel<-NULL
fullmodel1<-NULL
type3<-NULL
for (j in 1:length(y_type))
{if(y_type[j]==4 & is.null(w))
  {fullmodel1[[j]]<-coxph(y2[,j]~.,data=data.frame(cbind(x,pred)))
   fullmodel[[j]]<-summary(fullmodel1[[j]])}
 else if (y_type[j]==4)
  {fullmodel1[[j]]<-coxph(y2[,j]~.,data=data.frame(cbind(x,pred)), weights=w)
   fullmodel[[j]]<-summary(fullmodel1[[j]])}
 else
  {fullmodel1[[j]]<-glm(y2[,j]~.,data=data.frame(cbind(x,pred)),family=family1[[j]],weights=w) #use type three error to test the full model
   fullmodel[[j]]<-summary(fullmodel1[[j]])}
 type3[[j]]<-Anova(fullmodel1[[j]],type="III")
}
xname<-names(x)

xnames3<-rownames(type3[[1]])

P1<-matrix(NA,length(xnames3),ncol(y2))  #the type III for predictor and one mediator only model
rownames(P1)<-xnames3
colnames(P1)<-colnames(y2)
prednames<-colnames(pred)

covr.cont<-rep(F,length(contmed))
covr.bin<-rep(F,length(binmed))
covr.cat<-rep(F,length(catmed))

for (j in 1:ncol(y2))  #adjust for multivariate response
{if(testtype==2)
 {if(y_type[j]==4 & is.null(w))
   temp.fullmodel1<-coxph(y2[,j]~.,data=pred)
  else if (y_type[j]==4)
   temp.fullmodel1<-coxph(y2[,j]~.,data=pred,weights=w)
  else
   temp.fullmodel1<-glm(y2[,j]~.,data=pred,family=family1[[j]],weights=w) 
  temp.type3<-Anova(temp.fullmodel1,type="III")
  for (k in 1:ncol(pred))
   P1[grep(prednames[k],xnames3),j]<-temp.type3[rownames(temp.type3)==prednames[k],3]
 }

 if(!is.null(contmed))
  for (i in 1:length(contmed))
   if(testtype==1)
    covr.cont[i]<-ifelse(type3[[j]][xnames3==xname[contmed[i]],3]<alpha,T,covr.cont[i])
   else if(testtype==2)
   {temp.data<-cbind(x[,contmed[i]],pred)
    names(temp.data)<-c(xname[contmed[i]],names(pred))
    if(y_type[j]==4)
     temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
    else
     temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) 
    temp.type3<-Anova(temp.fullmodel1,type="III")
    temp.p1<-temp.type3[rownames(temp.type3)==xname[contmed[i]],3]
    P1[grep(xname[contmed[i]],xnames3),j]<-temp.p1
    covr.cont[i]<-ifelse(temp.p1<alpha,T,covr.cont[i])
   }

 if(!is.null(binmed))
  for (i in 1:length(binmed))
   if(testtype==1)
    covr.bin[i]<-ifelse(type3[[j]][xnames3==xname[binmed[i]],3]<alpha,T,covr.bin[i])
   else if(testtype==2)
   {temp.data<-cbind(x[,binmed[i]],pred)
    names(temp.data)<-c(xname[binmed[i]],names(pred))
    if(y_type[j]==4)
     temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
    else
     temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) 
   temp.type3<-Anova(temp.fullmodel1,type="III")
   temp.p1<-temp.type3[rownames(temp.type3)==xname[binmed[i]],3]
   P1[grep(xname[binmed[i]],xnames3),j]<-temp.p1
   covr.bin[i]<-ifelse(temp.p1<alpha,T,covr.bin[i])
  }
 
 if(!is.null(catmed))
  for (i in 1:length(catmed))
   if(testtype==1)
     covr.cat[i]<-ifelse(type3[[j]][xnames3==xname[catmed[i]],3]<alpha,T,covr.cat[i]) 
   else if(testtype==2)
    {temp.data<-cbind(x[,catmed[i]],pred)
     names(temp.data)<-c(xname[catmed[i]],names(pred))
     if(y_type[j]==4)
      temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
     else
      temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) #use type three error to test the full model
    temp.type3<-Anova(temp.fullmodel1,type="III")
    temp.p1<-temp.type3[rownames(temp.type3)==xname[catmed[i]],3]
    P1[grep(xname[catmed[i]],xnames3),j]<-temp.p1
    covr.cat[i]<-ifelse(temp.p1<alpha,T,covr.cat[i])
   }
} 

if(!is.null(contmed))
 {covr.cont<-ifelse(covr.cont|cont1,T,F)
  cont2<-cont1[covr.cont]
  contmed1<-contmed[covr.cont]}
if(!is.null(binmed))
 {covr.bin<-ifelse(covr.bin+bin1>0,T,F) 
  bin2<-bin1[covr.bin]
  binmed1<-binmed[covr.bin]}
if(!is.null(catmed))
 {covr.cat<-ifelse(covr.cat+cat1>0,T,F)
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
}
else {newx1<-x[,-cutx]
      if(sum(covr.cont)==0)
       contm1<-NULL  
      else 
       contm1<-colnum(contmed[covr.cont],cutx)
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
   {tempmodel<-summary(glm(newx1[,contm1[i]]~.,weights=w,data=pred)) #allowing multivariate predictors
    med.cont[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[contm1[i]])
    rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
   }
  med.cont<-ifelse(med.cont+cont2>0,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
    {tempmodel<-summary(glm(newx1[,binm1[i]]~.,weights=w,family="binomial",data=pred)) #allowing multivariate predictors
     med.bin[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
     rela_var<-c(rela_var,name_newx[binm1[i]])
     rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
    }
  med.bin<-ifelse(med.bin+bin2>0,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))  
   {temp.p<-NULL                                 #allowing multivariate predictors
    for (j in 1:ncol(pred))          
     temp.p<-c(temp.p,min(summary(glm(pred[,j]~newx1[,catm1[i]],weights=w,family="binomial"))$coef[-1,4]))
    med.cat[i]<-ifelse(min(temp.p)<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[catm1[i]])
    rela_p<-rbind(rela_p,temp.p)
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
  {tempmodel<-summary(glm(newx1[,contm1[i]]~.,weights=w,data=pred))
   med.cont[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
  }
 med.cont<-ifelse(med.cont|cont2,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
  {tempmodel<-summary(glm(newx1[,binm1[i]]~.,weights=w,family="binomial",data=pred))
   med.bin[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[binm1[i]])
   rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
  }    
  med.bin<-ifelse(med.bin|bin2,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))   
  {temp.p<-NULL
   for(j in 1:ncol(pred))
    temp.p<-c(temp.p,min(summary(glm(pred[,j]~newx1[,catm1[i]],weights=w))$coef[-1,4]))
   med.cat[i]<-ifelse(min(temp.p)<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[catm1[i]])
   rela_p<-rbind(rela_p,temp.p)
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
 rownames(rela_p)<-rela_var
 results<-list(x=newx2, dirx=pred, contm=contm2, catm=c(binm2,catm2),jointm=jointm,refy=refy,
              y=y2,y_type=y_type,fullmodel=fullmodel1,rela=rela_p,binpred=binpred,family1=family1,
              testtype=testtype,P1=P1,w=w)
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
 rownames(rela_p)<-rela_var
 results<-list(x=newx2,dirx=pred,contm=contm2,binm=binm2,catm=catm, jointm=jointm, refy=refy, y=y2, y_type=y_type,
               fullmodel=fullmodel1,rela=rela_p,binpred=binpred,family1=family1,
               testtype=testtype,P1=P1,w=w)
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
 covariates<-var.name[-c(object$contm,object$binm,t.catm)]
 tests<-NULL
 if(object$testtype==1)
  {for (j in 1:length(object$fullmodel))
    tests<-cbind(tests,Anova(object$fullmodel[[j]],type="III")[,3])
   temp<-rownames(Anova(object$fullmodel[[1]],type="III"))
   #tests<-cbind(tests,rep(NA,nrow(tests)))
 }
 else
  {tests<-object$P1 
   temp<-rownames(tests)
   #tests<-cbind(tests,rep(NA,nrow(tests)))
   }
  rownames(tests)<-temp
  temp.name<-rownames(object$rela)
  temp2<-matrix(NA,length(temp),ncol(object$rela))
 for (i in 1:nrow(object$rela))
   temp2[temp==temp.name[i],]<-object$rela[i,]
 tests<-cbind(tests,temp2)
 dimnames(tests)[[2]]<-c(paste("P-Value 1",colnames(object$y),sep="."), paste("P-Value 2", colnames(object$dirx),sep="."))
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
 dimnames(x$tests)[[2]]<-c(paste("P-Value 1",colnames(x$results$y),sep="."), paste("P-Value 2", colnames(x$results$dirx),sep="."))
 print(round(x$tests,3))
 cat("----\n *:mediator,-:joint mediator\n P-Value 1:Type-3 tests in the full model\n P-Value 2:Tests of relationship with the Predictor\n")
}

med<-function(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm,contm = data$contm, 
              catm = data$catm, jointm = data$jointm, allm = c(contm, catm), margin=1,
              n=20,seed=sample(1:1000,1), nonlinear=F, df=1, nu=0.001,D=3,distn=NULL,
              family1=data$family1,refy=rep(0,ncol(y)),binpred=data$binpred,x.new=x,
              pred.new=dirx,type=NULL, w=NULL, w.new=NULL)
{#for binary predictor
 med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                    catm = data$catm, jointm = data$jointm, allm = c(contm, catm), 
                    n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                    D=3,distn=NULL,family1=data$family1, #
                    biny=rep(F,ncol(y)),refy=rep(0,ncol(y)),surv=rep(F,ncol(y)),type=NULL, w=NULL) #
{if (is.null(allm))
  stop("Error: no potential mediator is specified")
  xnames<-colnames(x)
  pred_names<-colnames(dirx)  #
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  
  te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
  {te<-NULL
   for(m in 1:length(full.model))
    if(surv[m] & !is.null(best.iter1[m]))
     te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
    else if (surv[m])
     te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
    else
     te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
   te
  }
  
  med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
  {n3<-nrow(nom1)+nrow(nom0)
  marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
  new1<-nom1
  new1[,med]<-marg.m[1:nrow(nom1)]
  new0<-nom0
  new0[,med]<-marg.m[(nrow(nom1)+1):n3]
  dir.nom<-NULL
  for(m in 1:length(full.model))
   if(surv[m] & !is.null(best.iter1[m]))
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
   else if(surv[m])
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
   else
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
  dir.nom
  }
  
  med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand)  
  {if (length(med)==1)                       #added for the new program, when there is only one mediator
   {if(is.factor(nom1[,med]))              #added to control for one factor mediator
     marg.m<-as.factor(as.character(c(nom1[,med],nom0[,med])[temp.rand]))
    else
     marg.m<-c(nom1[,med],nom0[,med])[temp.rand]
   }        
  else                                         #added for the new program
   marg.m<-rbind(nom1[,med],nom0[,med])[temp.rand,]
  new1<-nom1
  new0<-nom0
  if(length(med)==1)                                       #added for the new program, when there is only one mediator
  {new1[,med]<-marg.m[1:nrow(new1)]                     #added for the new program 
   new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0))]}  #added for the new program
  else                                                     #added for the new program
  {new1[,med]<-marg.m[1:nrow(new1),]
   new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0)),]}
  dir.nom<-NULL
  for (m in 1:length(full.model))
   if(surv[m] & !is.null(best.iter1[m]))
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
   else if(surv[m])
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
   else
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
  dir.nom
  }
  
  med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
  {n3<-nrow(nom1)+nrow(nom0)
   temp.rand<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
   marg.m1<-temp.rand[1:nrow(nom1)]
   marg.m2<-temp.rand[(nrow(nom1)+1):n3]
   dir.nom<-rep(0,length(full.model))
   for (m in 1:length(full.model))
    for (i in levels(x[,med]))
    {new1<-nom1
     new1[1:dim(new1)[1],med]<-i
     new0<-nom0
     new0[1:dim(new0)[1],med]<-i
     p<-mean(temp.rand==i,na.rm=T)
     if(surv[m] & !is.null(best.iter1[m]))
       dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
     else if(surv[m])
       dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T))
     else
       dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
    }
  dir.nom
  }
  
  #1.fit the model
  x2<-cbind(x,dirx)
  colnames(x2)<-c(xnames,pred_names)
  full.model<-NULL
  best.iter1<-NULL
  if (nonlinear)
    for (j in 1:ncol(y)){
     if(biny[j])                     #recode y if y is binary
       y[,j]<-ifelse(y[,j]==refy[j],0,1)
     x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
     y1<-y[!is.na(y[,j]),j]
     full.model[[j]]<-suppressWarnings(gbm.fit(x1,y1, n.trees=200, interaction.depth=D, shrinkage=nu, w=w,
                                      distribution=distn[j],train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
     best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))
     while(full.model[[j]]$n.trees-best.iter1[j]<30){
      full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
      best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}}
  else
   for (j in 1:ncol(y)){
    if(biny[j])                     #recode y if y is binary
      y[,j]<-ifelse(y[,j]==refy[j],0,1)
    x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
    y1<-y[!is.na(y[,j]),j]
    if(surv[j])
     full.model[[j]]<-coxph(y1~., data=x1, weights=w)
    else
     full.model[[j]]<-glm(y1~., data=x1, family=family1[[j]], weights=w)
   }

  #2. prepare for the store of results
  set.seed(seed)
  te<-matrix(0,n,ncol(y)*ncol(dirx))
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  if(!is.null(jointm))
  {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
   dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
  }
  else
  {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))))
   dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)]),each=ncol(y)),sep=".")
  }
  denm<-rep(list(denm),ncol(dirx))
  ie<-denm
  #3. repeat to get the mediation effect
  for (k in 1:n)
  {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
    x0.temp<-apply(dirx==1,1,sum)==0  #indicator of the reference group
    x0<-x2[x0.temp,]
    if(is.null(w))
     {w1<-NULL
      w0<-NULL}
    else
      w0<-w[x0.temp]
    for (l in 1:ncol(dirx))  #l indicate the lth predictor
    {x1.2<-x2[dirx[,l]==1,]
     if(!is.null(w))
      w1<-w[dirx[,l]==1]
     #n3<-dim(x)[1] use the original size
     new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
     new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
     te[k,((l-1)*ncol(y)+1):(l*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
     temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
     denm[[l]][k,1:ncol(y)]<-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand) #add temp.rand
    j<-2
    #3.2 mediation effect from the continuous mediator
    if (!is.null(contm))
      for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type)
       j<-j+1}
    #3.3.mediation effect from the categorical mediator
    if (!is.null(catm))
      for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
      {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type)
       j<-j+1}
    #3.4 mediation effect from the joint mediators
    if (!is.null(jointm))
      for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
      {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
       denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand)
       j<-j+1}
    #3.5 get the indirect effects
    ie[[l]][k,]<-te[k,((l-1)*ncol(y)+1):(l*ncol(y))]-denm[[l]][k,]
    if(!is.null(jointm))
      dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
    else
      dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",names(x)[c(contm,catm)])
    }
  }
  names(denm)<-pred_names
  names(ie)<-pred_names
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),data=data)
  class(a)<-"med"
  return(a)
}

#for continous predictor
 med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                    catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                    nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=data$family1,
                    biny=(data$y_type==2),refy=rep(NA,ncol(y)),x.new=x,pred.new=dirx, surv=(data$y_type==4),type=NULL,w=NULL, w.new=NULL)
{if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  ynames<-colnames(y)
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
  
  anymissing<-function(vec) #return T if there is any missing in the vec
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
  col_mean<-function(col,n.row,w=NULL)
  {temp<-matrix(col,n.row)
  if(is.null(w))
    return(apply(temp,1,mean,na.rm=T))
  else
    return(apply(temp,1,weighted.mean,na.rm=T,w=w))}
  
  
  dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df,w) #give the model and residual of m given x
  {models<-NULL
  res<-NULL
  if(!is.null(catm))
  {for (i in 2:(catm$n+1))
    binm<-c(binm,catm[[i]])}
  if(nonlinear)
    {z<-NULL
     for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
  else
    z<-dirx
  j<-1
  if(!is.null(binm))
  {for(i in binm)
  {models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=binomial(link = "logit"),weights=w)
   res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
   j<-j+1}
  }
  for (i in contm)
  { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=gaussian(link="identity"),weights=w)
    res<-cbind(res,models[[j]]$res)
    j<-j+1
  }
  list(models=models,varmat=var(res))
  }
  
  
  sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df)  #added nonlinear and df to sim.xm
  {mult.norm<-function(mu,vari,n) 
   {if (nrow(vari)!=ncol(vari)) 
     result<-c("Error: Variance matrix is not square")  
    else if (length(mu)!=nrow(vari)) 
     result<-c("Error: length mu is not right!")  
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
  if(nonlinear)
  {z<-NULL
   for(i in 1:ncol(dirx))
    z<-cbind(z,ns(dirx[,i],df=df))}
   else
    z<-dirx
  
  binm1<-binm
  if(!is.null(catm))
  {for (i in 2:(catm$n+1))
    binm1<-c(binm1,catm[[i]])}
  
  if(!is.null(binm1))
    for (i in 1:length(binm1))
      means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z)))
  if(!is.null(contm))
    for (i in (length(binm1)+1):length(c(binm1,contm)))
      means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z)))
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
  
  if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
     tempm[[1]]<-NULL}
  else  tempm<-NULL
  allm<-unique(c(contm,binm,unlist(tempm)))
  
  nonmissing<-apply(cbind(y,x[,listm$single],dirx),1,anymissing)
  x<-x[nonmissing,]
  y<-data.frame(y[nonmissing,])
  colnames(y)<-ynames
  pred<-data.frame(dirx[nonmissing,])
  colnames(pred)<-pred_names
  w<-w[nonmissing]
  nonmissing1<-apply(cbind(x.new[,listm$single],pred.new),1,anymissing)
  x.new<-x.new[nonmissing1,]
  w.new<-w.new[nonmissing1]
  pred.new<-data.frame(pred.new[nonmissing1,])
  colnames(pred.new)<-pred_names
  
  #1.fit the model
  x2<-cbind(x,pred)
  colnames(x2)<-c(xnames,pred_names)
  full.model<-NULL
  best.iter1<-NULL
  for(j in 1:ncol(y)){
    if(biny[j])                     #recode y if y is binary
     y[,j]<-ifelse(y[,j]==refy[j],0,1)

  if(nonlinear)
  {full.model[[j]]<-suppressWarnings(gbm.fit(x2,y[,j], n.trees=200, interaction.depth=D, shrinkage=nu,w=w,
                                        distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
   best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))         
   while(full.model[[j]]$n.trees-best.iter1[j]<30){
    full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
    best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}
  }
  else
  {if(surv[j])
    full.model[[j]]<-coxph(y[,j]~., data=x2, weights=w)
   else
    full.model[[j]]<-glm(y[,j]~., data=x2, family=family1[[j]], weights=w)
  }
  }
  
  #2. prepare for the store of results
  set.seed(seed)
  n.new<-nrow(x.new)
#  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
  
  #3. get the joint distribution of m given x
  distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df,w)
  te1.0<-NULL
  denm1.0<-NULL

  n1<-dim(x)[1]

  #4. repeat to get the mediation effect
  for (l in 1:ncol(pred)) {
   denm1<-NULL
   te1<-NULL
   for (k in 1:n)
   {new0<-sim.xm(distmgivenx,x.new,pred.new,binm,contm,catm,nonlinear,df) #draw ms conditional on x.new
    temp.pred<-pred.new
    temp.pred[,l]<-temp.pred[,l]+margin
    new1<-sim.xm(distmgivenx,x.new,temp.pred,binm,contm,catm,nonlinear,df)  #draw from the conditional distribution of m given x
    new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
    new0<-cbind(new0,pred.new) 
    denm2<-NULL
    
    sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new)   #random sample from the original data
    
    temp.new1<-new1
    temp.new1[,allm]<-x.new[sample.temp[1:n.new],allm]
    temp.new0<-new0
    temp.new0[,allm]<-x.new[sample.temp[(n.new+1):(2*n.new)],allm]
    #4.0 get the direct effect
    for (m in 1:ncol(y))
     if(surv[m] & !is.null(best.iter1[m]))
       denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
     else if(surv[m])
       denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],temp.new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
     else
       denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)

    #4.1 get the te
    te0<-NULL
    for(m in 1:ncol(y))
     if(surv[m] & !is.null(best.iter1[m]))
      te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
     else if(surv[m])
      te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
     else
      te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
    te1<-cbind(te1,te0)
  
    #4.2 mediation effect from the single mediator
    if (!is.null(listm$single))
     for (i in 1:length(listm$single))
     {new1.nm<-new1
      new0.nm<-new0
      temp.m<-x.new[sample.temp,listm$single[i]]
      new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
      new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
      for(m in 1:ncol(y))
       if(surv[m] & !is.null(best.iter1[m]))
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
     }
  
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
    for (i in 2:(listm$multi[[1]]+1))
    {new1.nm<-new1
     new0.nm<-new0
     new1.nm[,listm$multi[[i]]]<-x.new[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
     new0.nm[,listm$multi[[i]]]<-x.new[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
     for(m in 1:col(y))
       if(surv[m] & !is.null(best.iter1[m]))
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
       else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
       else
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
    }
   denm1<-rbind(denm1,denm2)
  }
 denm1.0[[l]]<-denm1 
 te1.0[[l]]<-te1
 } 

   #4.4 get the indirect effects
  denm<-NULL
  te<-NULL
  ie<-NULL
  for (l in 1:ncol(pred))
   {denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
    te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
    te<-cbind(te,te0)
    temp1<-ncol(denm[[l]])/ncol(te0)
    temp2<-NULL
    for(temp in 1:temp1)
      temp2<-cbind(temp2,te0)
    ie[[l]]<-temp2-denm[[l]]
    if(!is.null(listm$multi)) 
     colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
     colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single]),each=ncol(y)),sep=".")
    if(!is.null(listm$multi))
     colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
     colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single]),each=ncol(y)),sep=".")
   }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  names(denm)<-pred_names
  names(ie)<-pred_names
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),pred.new=pred.new,w.new=w.new,data=data)
  class(a)<-"med"
  return(a)
}

 if(is.null(data)){
   surv=rep(F,ncol(y))
   biny=rep(F,ncol(y))
   if(is.null(distn))
     distn<-rep(NA,ncol(y))
   for(j in 1:ncol(y)) {
     if(class(y[,j])=="Surv"){
       surv[j]=T
       if(is.na(distn[j]))
         distn[j]="coxph"
       if(is.null(type) & nonlinear)
         type="response"
       else if (is.null(type))
         type="risk"
     }
     else if(is.character(y[,j]) | is.factor(y[,j]) | nlevels(as.factor(y[,j]))==2)
     {biny[j]=T
     if(is.na(family1[[j]]))
       family1[[j]] = binomial("logit")
     if(is.na(distn[j]))
       distn[j]="bernoulli" 
     if(!is.na(refy[j]))
       y[,j]<-ifelse(y[,j]==refy[j],0,1)
     else
       y[,j]<-ifelse(as.factor(y[,j])==levels(as.factor(y[,j]))[1],0,1)
     }
     else { 
       if(is.na(family1[[j]]))
         family1[[j]] = gaussian(link = "identity")
       if(is.na(distn[j]))
         distn[j]="gaussian" 
     }
   }
 }
 else
 {biny=data$y_type==2
 surv=data$y_type==4
 if(sum(surv)>0 & is.null(type) & nonlinear)
   type="response"
 else if (sum(surv)>0 & is.null(type))
   type="risk"
 if(is.null(distn))
   distn<-rep(NA,ncol(y))
 distn[is.na(distn) & data$y_type==2]="bernoulli"
 distn[is.na(distn) & data$y_type==4]="coxph"
 distn[is.na(distn) & data$y_type==1]="gaussian"
 }
 
  if(binpred)
    a<-med.binx(data=data, x=x, y=y, dirx=dirx, contm = contm, 
                catm=catm, jointm=jointm, allm=allm, n=n,seed=seed,
                nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,
                biny=biny,refy=refy,surv=surv,type=type,w)
  else
    a<-med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                 catm=catm, jointm=jointm, margin=margin, n=n, seed=seed, 
                 nonlinear=nonlinear, df=df, nu=nu,D=D, distn=distn, 
                 family1=family1,biny=biny,refy=refy,x.new=x.new,pred.new=pred.new,surv=surv,type=type,w,w.new)
  return(a)
}
  

print.med<-function(x,...,digit=4)
{for(l in 1:length(x$ie))
 {cat("\n\nFor the predictor",names(x$ie)[l],":\n")
  cat(" The estimated total effect:")
  if(is.null(x$w.new))
    print(mean(x$te[,l],na.rm=T),digit)
  else
    print(round(weighted.mean(x$te[,l],na.rm=T,w=x$w.new),digit))
  cat("\n The estimated indirect effect:\n")
  if(is.null(x$w.new))
     print(round(apply(x$ie[[l]],2,mean,na.rm=T),digit))
  else
     print(round(apply(x$ie[[l]],2,weighted.mean,na.rm=T,w=x$w.new),digit))}
}


boot.med<-function(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,catm=data$catm,
                        jointm=data$jointm,margin=1,n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,
                        D=3,distn=NULL,family1=data$family1,n2=50,w=rep(1,nrow(x)),
                        refy=NULL,x.new=x,pred.new=dirx,binpred=data$binpred,type=NULL,w.new=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),surv=(data$y_type==4),type)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, allm = c(contm, catm), 
                     n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn=NULL,family1=data$family1, #
                     biny=rep(F,ncol(y)),refy=rep(0,ncol(y)),surv=rep(F,ncol(y)),type=NULL, w=NULL) #
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    pred_names<-colnames(dirx)  #
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {te<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if (surv[m])
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      te
    }
    
    med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {if(nrow(nom1)==0 | nrow(nom0)==0)            #in case a group of predictors has very few sample size
      return(rep(NA,length(full.model)))
    n3<-nrow(nom1)+nrow(nom0)
    marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    new1<-nom1
    new1[,med]<-marg.m[1:nrow(nom1)]
    new0<-nom0
    new0[,med]<-marg.m[(nrow(nom1)+1):n3]
    dir.nom<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
    else if(surv[m])
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand)  
    {if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-as.factor(as.character(c(nom1[,med],nom0[,med])[temp.rand]))
    else
      marg.m<-c(nom1[,med],nom0[,med])[temp.rand]
    }        
      else                                         #added for the new program
        marg.m<-rbind(nom1[,med],nom0[,med])[temp.rand,]
      new1<-nom1
      new0<-nom0
      if(length(med)==1)                                       #added for the new program, when there is only one mediator
      {new1[,med]<-marg.m[1:nrow(new1)]                     #added for the new program 
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0))]}  #added for the new program
      else                                                     #added for the new program
      {new1[,med]<-marg.m[1:nrow(new1),]
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0)),]}
      dir.nom<-NULL
      for (m in 1:length(full.model))
        if(surv[m] & !is.null(best.iter1[m]))
          dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if(surv[m])
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      dir.nom
    }
    
    med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    temp.rand<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    marg.m1<-temp.rand[1:nrow(nom1)]
    marg.m2<-temp.rand[(nrow(nom1)+1):n3]
    dir.nom<-rep(0,length(full.model))
    for (m in 1:length(full.model))
      for (i in levels(x[,med]))
      {new1<-nom1
      new1[1:dim(new1)[1],med]<-i
      new0<-nom0
      new0[1:dim(new0)[1],med]<-i
      p<-mean(temp.rand==i,na.rm=T)
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
      else if(surv[m])
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T))
      else
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
      }
    dir.nom
    }
    
    #1.fit the model
    x2<-cbind(x,dirx)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    if (nonlinear)
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        full.model[[j]]<-suppressWarnings(gbm.fit(x1,y1, n.trees=200, interaction.depth=D, shrinkage=nu, w=w,
                                                  distribution=distn[j],train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))
        while(full.model[[j]]$n.trees-best.iter1[j]<30){
          full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
          best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}}
    else
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        if(surv[j])
          full.model[[j]]<-coxph(y1~., data=x1, weights=w)
        else
          full.model[[j]]<-glm(y1~., data=x1, family=family1[[j]], weights=w)
      }
    
    #2. prepare for the store of results
    set.seed(seed)
    te<-matrix(0,n,ncol(y)*ncol(dirx))
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    if(!is.null(jointm))
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
    }
    else
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)]),each=ncol(y)),sep=".")
    }
    denm<-rep(list(denm),ncol(dirx))
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x0.temp<-apply(dirx==1,1,sum)==0  #indicator of the reference group
      x0<-x2[x0.temp,]
      if(is.null(w))
      {w1<-NULL
      w0<-NULL}
      else
        w0<-w[x0.temp]
      for (l in 1:ncol(dirx))  #l indicate the lth predictor
      {x1.2<-x2[dirx[,l]==1,]
      if(!is.null(w))
        w1<-w[dirx[,l]==1]
      #n3<-dim(x)[1] use the original size
      new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
      new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
      te[k,((l-1)*ncol(y)+1):(l*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
      temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
      denm[[l]][k,1:ncol(y)]<-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand) #add temp.rand
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
        denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand)
        j<-j+1}
      #3.5 get the indirect effects
      ie[[l]][k,]<-te[k,((l-1)*ncol(y)+1):(l*ncol(y))]-denm[[l]][k,]
      if(!is.null(jointm))
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",names(x)[c(contm,catm)])
      }
    }
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),data=data)
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  if(is.null(jointm))
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))))
   ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))))
   dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
   colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
   rownames(ie1)<-pred_names}
  else 
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))+jointm[[1]]))
   dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
   ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))+jointm[[1]]))
   dimnames(ie1)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
   rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  temp<-med.binx(data=NULL,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type,w=w)
  te[1,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
   ie1[l,]<-apply(temp$ie[[l]],2,mean)}  #first row is the estimated value
  de[1,]<-apply(temp.1,2,mean,na.rm=T)
  model<-temp$model
  best.iter<-temp$best.iter
  
  for (i in 1:n2)
  {boots<-sample(1:nrow(x),replace=T,prob=w)
   x1<-x[boots,]
   y1<-data.frame(y[boots,])
   pred1<-data.frame(dirx[boots,])
   temp<-med.binx(data=NULL,x=x1, y=y1, dirx=pred1, contm=contm, catm=catm,jointm=jointm,allm=allm,n=n,seed=seed+i,
                  nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,biny=biny,refy=refy,surv=surv,type=type,w=NULL)
   te[1+i,]<-apply(temp$te,2,mean,na.rm=T)
   temp.1<-NULL
   for (l in 1:nx)
   {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
    ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)}  #first row is the estimated value
   de[1+i,]<-apply(temp.1,2,mean,na.rm=T)
   print(i)
  }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")

  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),
                         x.new=x,pred.new=dirx,surv,type,w.new=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=data$family1,
                      biny=(data$y_type==2),refy=rep(NA,ncol(y)),x.new=x,pred.new=dirx, surv=(data$y_type==4),type=NULL,w=NULL, w.new=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    pred_names<-colnames(dirx)
    ynames<-colnames(y)
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
    
    anymissing<-function(vec) #return T if there is any missing in the vec
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row,w=NULL)
    {temp<-matrix(col,n.row)
    if(is.null(w))
      return(apply(temp,1,mean,na.rm=T))
    else
      return(apply(temp,1,weighted.mean,na.rm=T,w=w))}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df,w) #give the model and residual of m given x
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=binomial(link = "logit"),weights=w)
      res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
      j<-j+1}
    }
    for (i in contm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=gaussian(link="identity"),weights=w)
      res<-cbind(res,models[[j]]$res)
      j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df)  #added nonlinear and df to sim.xm
    {mult.norm<-function(mu,vari,n) 
    {if (nrow(vari)!=ncol(vari)) 
      result<-c("Error: Variance matrix is not square")  
    else if (length(mu)!=nrow(vari)) 
      result<-c("Error: length mu is not right!")  
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
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    
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
    
    if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
    tempm[[1]]<-NULL}
    else  tempm<-NULL
    allm<-unique(c(contm,binm,unlist(tempm)))
    
    nonmissing<-apply(cbind(y,x[,listm$single],dirx),1,anymissing)
    x<-x[nonmissing,]
    y<-data.frame(y[nonmissing,])
    colnames(y)<-ynames
    pred<-data.frame(dirx[nonmissing,])
    colnames(pred)<-pred_names
    w<-w[nonmissing]
    nonmissing1<-apply(cbind(x.new[,listm$single],pred.new),1,anymissing)
    x.new<-x.new[nonmissing1,]
    w.new<-w.new[nonmissing1]
    pred.new<-data.frame(pred.new[nonmissing1,])
    colnames(pred.new)<-pred_names
    
    #1.fit the model
    x2<-cbind(x,pred)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    for(j in 1:ncol(y)){
      if(biny[j])                     #recode y if y is binary
        y[,j]<-ifelse(y[,j]==refy[j],0,1)
      
      if(nonlinear)
      {full.model[[j]]<-suppressWarnings(gbm.fit(x2,y[,j], n.trees=200, interaction.depth=D, shrinkage=nu,w=w,
                                                 distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
      best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))         
      while(full.model[[j]]$n.trees-best.iter1[j]<30){
        full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}
      }
      else
      {if(surv[j])
        full.model[[j]]<-coxph(y[,j]~., data=x2, weights=w)
      else
        full.model[[j]]<-glm(y[,j]~., data=x2, family=family1[[j]], weights=w)
      }
    }
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    #  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df,w)
    te1.0<-NULL
    denm1.0<-NULL
    
    n1<-dim(x)[1]
    
    #4. repeat to get the mediation effect
    for (l in 1:ncol(pred)) {
      denm1<-NULL
      te1<-NULL
      for (k in 1:n)
      {new0<-sim.xm(distmgivenx,x.new,pred.new,binm,contm,catm,nonlinear,df) #draw ms conditional on x.new
      temp.pred<-pred.new
      temp.pred[,l]<-temp.pred[,l]+margin
      new1<-sim.xm(distmgivenx,x.new,temp.pred,binm,contm,catm,nonlinear,df)  #draw from the conditional distribution of m given x
      new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
      new0<-cbind(new0,pred.new) 
      denm2<-NULL
      
      sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new)   #random sample from the original data
      
      temp.new1<-new1
      temp.new1[,allm]<-x.new[sample.temp[1:n.new],allm]
      temp.new0<-new0
      temp.new0[,allm]<-x.new[sample.temp[(n.new+1):(2*n.new)],allm]
      #4.0 get the direct effect
      for (m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],temp.new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)
      
      #4.1 get the te
      te0<-NULL
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
      te1<-cbind(te1,te0)
      
      #4.2 mediation effect from the single mediator
      if (!is.null(listm$single))
        for (i in 1:length(listm$single))
        {new1.nm<-new1
        new0.nm<-new0
        temp.m<-x.new[sample.temp,listm$single[i]]
        new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
        new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
        for(m in 1:ncol(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      
      #4.3.mediation effect from the joint mediator
      if (!is.null(listm$multi))
        for (i in 2:(listm$multi[[1]]+1))
        {new1.nm<-new1
        new0.nm<-new0
        new1.nm[,listm$multi[[i]]]<-x.new[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
        new0.nm[,listm$multi[[i]]]<-x.new[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
        for(m in 1:col(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      denm1<-rbind(denm1,denm2)
      }
      denm1.0[[l]]<-denm1 
      te1.0[[l]]<-te1
    } 
    
    #4.4 get the indirect effects
    denm<-NULL
    te<-NULL
    ie<-NULL
    for (l in 1:ncol(pred))
    {denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
    te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
    te<-cbind(te,te0)
    temp1<-ncol(denm[[l]])/ncol(te0)
    temp2<-NULL
    for(temp in 1:temp1)
      temp2<-cbind(temp2,te0)
    ie[[l]]<-temp2-denm[[l]]
    if(!is.null(listm$multi)) 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single]),each=ncol(y)),sep=".")
    if(!is.null(listm$multi))
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single]),each=ncol(y)),sep=".")
    }
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),pred.new=pred.new,w.new=w.new,data=data)
    class(a)<-"med"
    return(a)
  }
  
  anymissing<-function(vec)
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")

xnames<-colnames(x)
pred_names<-colnames(dirx)
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

ny=ncol(y)
nx=ncol(dirx)
te<-matrix(0,n2+1,ny*nx)
de<-matrix(0,n2+1,ny*nx)
mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
ie<-matrix(0,n2,ny*(1+length(listm$single)+mul))   #added in the new program
ie1<-matrix(0,nx,ny*(1+length(listm$single)+mul))   #added in the new program
if(!is.null(listm$multi))
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
   colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
   rownames(ie1)<-pred_names}
else 
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
   colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
   rownames(ie1)<-pred_names}
ie<-rep(list(ie),nx)
names(ie)<-pred_names

temp<-med.contx(data=NULL,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                refy=refy,x.new=x.new,pred.new=pred.new, surv,type,w=w,w.new=w.new)

temp.1<-NULL
for (l in 1:nx)
 temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
if(is.null(w.new))
{te[1,]<-apply(temp$te,2,mean,na.rm=T)
 de[1,]<-apply(temp.1,2,mean,na.rm=T) 
 for (l in 1:nx)
   ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
}
else
{te[1,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
 de[1,]<-apply(temp$denm[,1],2,weighted.mean,na.rm=T,w=w.new) 
 for (l in 1:nx)
   ie1[l,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new)  #first row is the estimated value
}

te1<-NULL                      #to store the mediation effects on predictor
de1<-NULL
ie2<-rep(list(NULL),nx)
names(ie2)<-pred_names
model<-temp$model
for (i in 1:n2)
{boots<-sample(1:nrow(x),replace=T, prob=w)
 x1<-x[boots,]
 y1<-data.frame(y[boots,])
 dirx1<-data.frame(dirx[boots,])
 temp<-med.contx(data=NULL,x=x1,y=y1,dirx=dirx1,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                 margin=margin,n=n,seed=seed+i,nonlinear=nonlinear,df=df,nu=nu,D=D,
                 distn=distn,family1=family1,biny=biny,refy=refy,x.new=x.new,pred.new=pred.new,surv=surv,type=type) #added to the new codel, change the seed to make different results
 temp.1<-NULL
 for (l in 1:nx)
   temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
 if(is.null(w.new))
   {te[1+i,]<-apply(temp$te,2,mean,na.rm=T)
    de[1+i,]<-apply(temp.1,2,mean,na.rm=T)
    for (l in 1:nx)
      ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
   }
else
{te[1+i,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
 de[1+i,]<-apply(temp$denm[,1],weighted.mean,na.rm=T,w=w.new)
 for (l in 1:nx)
   ie[[l]][i,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T)  #first row is the estimated value
}
te1<-cbind(te1,temp$te)
de1<-cbind(de1,temp.1)
for (l in 1:nx)
  ie2[[l]]<-rbind(ie2[[l]],temp$ie[[l]])
print(i)
}
colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
missing.pred.new<-apply(pred.new,1,anymissing)
pred.new<-data.frame(pred.new[!missing.pred.new,])

a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model,
        data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
        boot.detail=list(pred.new=pred.new,te1=te1,de1=de1,ie1=ie2),w.new=w.new)
class(a)<-"mma"
return(a)
}


if(is.null(data)){
  surv=rep(F,ncol(y))
  biny=rep(F,ncol(y))
  if(is.null(distn))
    distn<-rep(NA,ncol(y))
  for(j in 1:ncol(y)) {
    if(class(y[,j])=="Surv"){
      surv[j]=T
      if(is.na(distn[j]))
        distn[j]="coxph"
      if(is.null(type) & nonlinear)
        type="response"
      else if (is.null(type))
        type="risk"
    }
    else if(is.character(y[,j]) | is.factor(y[,j]) | nlevels(as.factor(y[,j]))==2)
    {biny[j]=T
    if(is.na(family1[[j]]))
      family1[[j]] = binomial("logit")
    if(is.na(distn[j]))
      distn[j]="bernoulli" 
    if(!is.na(refy[j]))
      y[,j]<-ifelse(y[,j]==refy[j],0,1)
    else
      y[,j]<-ifelse(as.factor(y[,j])==levels(as.factor(y[,j]))[1],0,1)
    }
    else { 
      if(is.na(family1[[j]]))
        family1[[j]] = gaussian(link = "identity")
      if(is.na(distn[j]))
        distn[j]="gaussian" 
    }
  }
}
else
{biny=data$y_type==2
surv=data$y_type==4
if(sum(surv)>0 & is.null(type) & nonlinear)
  type="response"
else if (sum(surv)>0 & is.null(type))
  type="risk"
if(is.null(distn))
  distn<-rep(NA,ncol(y))
distn[is.na(distn) & data$y_type==2]="bernoulli"
distn[is.na(distn) & data$y_type==4]="coxph"
distn[is.na(distn) & data$y_type==1]="gaussian"
}

if(binpred)
  a<-boot.med.binx(data=data,x=x, y=y,dirx=dirx,contm=contm,catm=catm,
                             jointm=jointm,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,
                             D=D,distn=distn,family1=family1,
                             w=w,biny=biny,refy=rep(0,ncol(y)),surv,type)
else
  a<-boot.med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                    catm=catm, jointm=jointm, margin = margin, n = n, seed = seed, 
                    nonlinear = nonlinear, df = df, nu = nu, D = D, distn = distn, 
                    family1 = family1, n2 = n2,w=w,
                    biny=biny,refy=rep(0,ncol(y)),x.new=x.new,pred.new=pred.new, surv,type,w.new)

return(a)
}

    
  

mma<-function(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,
              catmed=NULL,catref=NULL,jointm=NULL,refy=rep(NA,ncol(data.frame(y))),
              predref=NULL,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
              nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=as.list(rep(NA,ncol(data.frame(y)))),
              n2=50,w=rep(1,nrow(x)), testtype=1, x.new=NULL, pred.new=NULL, type=NULL,w.new=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),surv=(data$y_type==4),type)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, allm = c(contm, catm), 
                     n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn=NULL,family1=data$family1, #
                     biny=rep(F,ncol(y)),refy=rep(0,ncol(y)),surv=rep(F,ncol(y)),type=NULL, w=NULL) #
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    pred_names<-colnames(dirx)  #
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {te<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if (surv[m])
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      te
    }
    
    med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {if(nrow(nom1)==0 | nrow(nom0)==0)            #in case a group of predictors has very few sample size
      return(rep(NA,length(full.model)))
      n3<-nrow(nom1)+nrow(nom0)
      marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
      new1<-nom1
      new1[,med]<-marg.m[1:nrow(nom1)]
      new0<-nom0
      new0[,med]<-marg.m[(nrow(nom1)+1):n3]
      dir.nom<-NULL
      for(m in 1:length(full.model))
        if(surv[m] & !is.null(best.iter1[m]))
          dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if(surv[m])
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      dir.nom
    }
    
    med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand)  
    {if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-as.factor(as.character(c(nom1[,med],nom0[,med])[temp.rand]))
    else
      marg.m<-c(nom1[,med],nom0[,med])[temp.rand]
    }        
      else                                         #added for the new program
        marg.m<-rbind(nom1[,med],nom0[,med])[temp.rand,]
      new1<-nom1
      new0<-nom0
      if(length(med)==1)                                       #added for the new program, when there is only one mediator
      {new1[,med]<-marg.m[1:nrow(new1)]                     #added for the new program 
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0))]}  #added for the new program
      else                                                     #added for the new program
      {new1[,med]<-marg.m[1:nrow(new1),]
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0)),]}
      dir.nom<-NULL
      for (m in 1:length(full.model))
        if(surv[m] & !is.null(best.iter1[m]))
          dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if(surv[m])
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      dir.nom
    }
    
    med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    temp.rand<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    marg.m1<-temp.rand[1:nrow(nom1)]
    marg.m2<-temp.rand[(nrow(nom1)+1):n3]
    dir.nom<-rep(0,length(full.model))
    for (m in 1:length(full.model))
      for (i in levels(x[,med]))
      {new1<-nom1
      new1[1:dim(new1)[1],med]<-i
      new0<-nom0
      new0[1:dim(new0)[1],med]<-i
      p<-mean(temp.rand==i,na.rm=T)
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
      else if(surv[m])
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T))
      else
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
      }
    dir.nom
    }
    
    #1.fit the model
    x2<-cbind(x,dirx)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    if (nonlinear)
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        full.model[[j]]<-suppressWarnings(gbm.fit(x1,y1, n.trees=200, interaction.depth=D, shrinkage=nu, w=w,
                                                  distribution=distn[j],train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))
        while(full.model[[j]]$n.trees-best.iter1[j]<30){
          full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
          best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}}
    else
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        if(surv[j])
          full.model[[j]]<-coxph(y1~., data=x1, weights=w)
        else
          full.model[[j]]<-glm(y1~., data=x1, family=family1[[j]], weights=w)
      }
    
    #2. prepare for the store of results
    set.seed(seed)
    te<-matrix(0,n,ncol(y)*ncol(dirx))
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    if(!is.null(jointm))
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
    }
    else
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)]),each=ncol(y)),sep=".")
    }
    denm<-rep(list(denm),ncol(dirx))
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x0.temp<-apply(dirx==1,1,sum)==0  #indicator of the reference group
      x0<-x2[x0.temp,]
      if(is.null(w))
      {w1<-NULL
      w0<-NULL}
      else
        w0<-w[x0.temp]
      for (l in 1:ncol(dirx))  #l indicate the lth predictor
      {x1.2<-x2[dirx[,l]==1,]
      if(!is.null(w))
        w1<-w[dirx[,l]==1]
      #n3<-dim(x)[1] use the original size
      new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
      new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
      te[k,((l-1)*ncol(y)+1):(l*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
      temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
      denm[[l]][k,1:ncol(y)]<-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand) #add temp.rand
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
        denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand)
        j<-j+1}
      #3.5 get the indirect effects
      ie[[l]][k,]<-te[k,((l-1)*ncol(y)+1):(l*ncol(y))]-denm[[l]][k,]
      if(!is.null(jointm))
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",names(x)[c(contm,catm)])
      }
    }
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),data=data)
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  if(is.null(jointm))
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))))
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie1)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  temp<-med.binx(data=NULL,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type,w=w)
  te[1,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie1[l,]<-apply(temp$ie[[l]],2,mean)}  #first row is the estimated value
  de[1,]<-apply(temp.1,2,mean,na.rm=T)
  model<-temp$model
  best.iter<-temp$best.iter
  
  for (i in 1:n2)
  {boots<-sample(1:nrow(x),replace=T,prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  pred1<-data.frame(dirx[boots,])
  temp<-med.binx(data=NULL,x=x1, y=y1, dirx=pred1, contm=contm, catm=catm,jointm=jointm,allm=allm,n=n,seed=seed+i,
                 nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,biny=biny,refy=refy,surv=surv,type=type,w=NULL)
  te[1+i,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)}  #first row is the estimated value
  de[1+i,]<-apply(temp.1,2,mean,na.rm=T)
  print(i)
  }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),
                         x.new=x,pred.new=dirx,surv,type,w.new=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=data$family1,
                      biny=(data$y_type==2),refy=rep(NA,ncol(y)),x.new=x,pred.new=dirx, surv=(data$y_type==4),type=NULL,w=NULL, w.new=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    pred_names<-colnames(dirx)
    ynames<-colnames(y)
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
    
    anymissing<-function(vec) #return T if there is any missing in the vec
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row,w=NULL)
    {temp<-matrix(col,n.row)
    if(is.null(w))
      return(apply(temp,1,mean,na.rm=T))
    else
      return(apply(temp,1,weighted.mean,na.rm=T,w=w))}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df,w) #give the model and residual of m given x
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=binomial(link = "logit"),weights=w)
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=gaussian(link="identity"),weights=w)
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df)  #added nonlinear and df to sim.xm
    {mult.norm<-function(mu,vari,n) 
    {if (nrow(vari)!=ncol(vari)) 
      result<-c("Error: Variance matrix is not square")  
    else if (length(mu)!=nrow(vari)) 
      result<-c("Error: length mu is not right!")  
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
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    
    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z)))
    
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z)))
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
    
    if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
    tempm[[1]]<-NULL}
    else  tempm<-NULL
    allm<-unique(c(contm,binm,unlist(tempm)))
    
    nonmissing<-apply(cbind(y,x[,listm$single],dirx),1,anymissing)
    x<-x[nonmissing,]
    y<-data.frame(y[nonmissing,])
    colnames(y)<-ynames
    pred<-data.frame(dirx[nonmissing,])
    colnames(pred)<-pred_names
    w<-w[nonmissing]
    nonmissing1<-apply(cbind(x.new[,listm$single],pred.new),1,anymissing)
    x.new<-x.new[nonmissing1,]
    w.new<-w.new[nonmissing1]
    pred.new<-data.frame(pred.new[nonmissing1,])
    colnames(pred.new)<-pred_names
    
    #1.fit the model
    x2<-cbind(x,pred)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    for(j in 1:ncol(y)){
      if(biny[j])                     #recode y if y is binary
        y[,j]<-ifelse(y[,j]==refy[j],0,1)
      
      if(nonlinear)
      {full.model[[j]]<-suppressWarnings(gbm.fit(x2,y[,j], n.trees=200, interaction.depth=D, shrinkage=nu,w=w,
                                                 distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
      best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))         
      while(full.model[[j]]$n.trees-best.iter1[j]<30){
        full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}
      }
      else
      {if(surv[j])
        full.model[[j]]<-coxph(y[,j]~., data=x2, weights=w)
      else
        full.model[[j]]<-glm(y[,j]~., data=x2, family=family1[[j]], weights=w)
      }
    }
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    #  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df,w)
    te1.0<-NULL
    denm1.0<-NULL
    
    n1<-dim(x)[1]
    
    #4. repeat to get the mediation effect
    for (l in 1:ncol(pred)) {
      denm1<-NULL
      te1<-NULL
      for (k in 1:n)
      {new0<-sim.xm(distmgivenx,x.new,pred.new,binm,contm,catm,nonlinear,df) #draw ms conditional on x.new
      temp.pred<-pred.new
      temp.pred[,l]<-temp.pred[,l]+margin
      new1<-sim.xm(distmgivenx,x.new,temp.pred,binm,contm,catm,nonlinear,df)  #draw from the conditional distribution of m given x
      new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
      new0<-cbind(new0,pred.new) 
      denm2<-NULL
      
      sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new)   #random sample from the original data
      
      temp.new1<-new1
      temp.new1[,allm]<-x.new[sample.temp[1:n.new],allm]
      temp.new0<-new0
      temp.new0[,allm]<-x.new[sample.temp[(n.new+1):(2*n.new)],allm]
      #4.0 get the direct effect
      for (m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],temp.new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)
      
      #4.1 get the te
      te0<-NULL
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
      te1<-cbind(te1,te0)
      
      #4.2 mediation effect from the single mediator
      if (!is.null(listm$single))
        for (i in 1:length(listm$single))
        {new1.nm<-new1
        new0.nm<-new0
        temp.m<-x.new[sample.temp,listm$single[i]]
        new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
        new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
        for(m in 1:ncol(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      
      #4.3.mediation effect from the joint mediator
      if (!is.null(listm$multi))
        for (i in 2:(listm$multi[[1]]+1))
        {new1.nm<-new1
        new0.nm<-new0
        new1.nm[,listm$multi[[i]]]<-x.new[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
        new0.nm[,listm$multi[[i]]]<-x.new[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
        for(m in 1:col(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      denm1<-rbind(denm1,denm2)
      }
      denm1.0[[l]]<-denm1 
      te1.0[[l]]<-te1
    } 
    
    #4.4 get the indirect effects
    denm<-NULL
    te<-NULL
    ie<-NULL
    for (l in 1:ncol(pred))
    {denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
    te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
    te<-cbind(te,te0)
    temp1<-ncol(denm[[l]])/ncol(te0)
    temp2<-NULL
    for(temp in 1:temp1)
      temp2<-cbind(temp2,te0)
    ie[[l]]<-temp2-denm[[l]]
    if(!is.null(listm$multi)) 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single]),each=ncol(y)),sep=".")
    if(!is.null(listm$multi))
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single]),each=ncol(y)),sep=".")
    }
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),x.new=x.new,w.new=w.new,data=data)
    class(a)<-"med"
    return(a)
  }
  
  anymissing<-function(vec)
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
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
  
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2,ny*(1+length(listm$single)+mul))   #added in the new program
  ie1<-matrix(0,nx,ny*(1+length(listm$single)+mul))   #added in the new program
  if(!is.null(listm$multi))
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  
  temp<-med.contx(data=NULL,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,pred.new=pred.new, surv,type,w=w,w.new=w.new)
  
  temp.1<-NULL
  for (l in 1:nx)
    temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  if(is.null(w.new))
  {te[1,]<-apply(temp$te,2,mean,na.rm=T)
  de[1,]<-apply(temp.1,2,mean,na.rm=T) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
  }
  else
  {te[1,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
  de[1,]<-apply(temp$denm[,1],2,weighted.mean,na.rm=T,w=w.new) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new)  #first row is the estimated value
  }
  
  te1<-NULL                      #to store the mediation effects on predictor
  de1<-NULL
  ie2<-rep(list(NULL),nx)
  names(ie2)<-pred_names
  model<-temp$model
  for (i in 1:n2)
  {boots<-sample(1:nrow(x),replace=T, prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  dirx1<-data.frame(dirx[boots,])
  temp<-med.contx(data=NULL,x1,y1,dirx1,binm,contm,catm,jointm, margin,n,seed+i,nonlinear,df,nu,D,
                  distn,family1,biny,refy,x.new,pred.new,surv,type) #added to the new codel, change the seed to make different results
  temp.1<-NULL
  for (l in 1:nx)
    temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  if(is.null(w.new))
  {te[1+i,]<-apply(temp$te,2,mean,na.rm=T)
  de[1+i,]<-apply(temp.1,2,mean,na.rm=T)
  for (l in 1:nx)
    ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
  }
  else
  {te[1+i,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
  de[1+i,]<-apply(temp$denm[,1],weighted.mean,na.rm=T,w=w.new)
  for (l in 1:nx)
    ie[[l]][i,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T)  #first row is the estimated value
  }
  te1<-cbind(te1,temp$te)
  de1<-cbind(de1,temp.1)
  for (l in 1:nx)
    ie2[[l]]<-rbind(ie2[[l]],temp$ie[[l]])
  print(i)
  }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  missing.pred.new<-apply(pred.new,1,anymissing)
  pred.new<-data.frame(pred.new[missing.pred.new,])
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(pred.new=pred.new,te1=te1,de1=de1,ie1=ie2),w.new=w.new)
  class(a)<-"mma"
  return(a)
}

data<-data.org(x=x,y=y,pred=pred,mediator=mediator,contmed=contmed,binmed=binmed,
               binref=binref,catmed=catmed,catref=catref,jointm=jointm,refy=refy,family1=family1,
               predref=predref,alpha=alpha,alpha2=alpha2,testtype=testtype, w=w)
biny=data$y_type==2
surv=data$y_type==4
if(sum(surv)>0 & is.null(type) & nonlinear)
  type="response"
else if (sum(surv)>0 & is.null(type))
  type="risk"
if(is.null(distn))
  distn<-rep(NA,ncol(data$y))
distn[is.na(distn) & data$y_type==2]="bernoulli"
distn[is.na(distn) & data$y_type==4]="coxph"
distn[is.na(distn) & data$y_type==1]="gaussian"

binpred<-data$binpred
family1<-data$family1

if(binpred) 
  {result<-boot.med.binx(data=data,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,
                         w=w,biny=biny,refy=rep(0,ncol(data$y)),surv=surv,type=type)
  }
 else
  {if(is.null(pred.new))
    result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                           D=D,distn=distn,family1=family1,n2=n2,w=w,biny=biny,refy=rep(0,ncol(data$y)),surv=surv,type=type)
   else
    result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                           D=D,distn=distn,family1=family1,n2=n2,w=w,biny=biny,refy=0, x.new=x.new, pred.new=pred.new,surv=surv,type=type,w.new=w.new)
  }
 result
}


#classes and methods for mma
print.mma<-function(x,...,digit=3)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$model$MART)
   cat ("MART\n")
 else cat("GLM\n")
 print(x$e,digit=digit)
}


summary.mma<-function(object,...,alpha=0.05,plot=TRUE,RE=FALSE,quant=T)
{x<-object
 ny<-ncol(x$data$y)
 nx<-ncol(x$data$dirx)
 temp1<-x$boots
 temp2<-x$est
 
 temp3<-x$boots   #calculate the RE
 temp3$de<-temp3$de/temp3$te
 nie<-ncol(temp3$ie[[1]])/ny
 k<-0
 for(l in 1:nx)
   {temp.te<-as.matrix(temp3$te)[,(k+1):(k+ny)]
    temp.te1<-do.call(cbind, replicate(nie,temp.te,simplify=F))
    temp3$ie[[l]]<-temp3$ie[[l]]/temp.te1
    k<-k+ny
   }
 temp4<-x$est
 temp.te<-matrix(temp4$te,nx,ny,byrow=T)
 temp.te1<-do.call(cbind, replicate(nie,temp.te,simplify=F))
 temp4$ie<-temp4$ie/temp.te1
 temp4$de<-temp4$de/temp4$te
 
 a1<-alpha/2
 a2<-1-a1
 b1<-qnorm(a1)
 b2<-qnorm(a2)
 ie<-NULL
 for (l in 1:nx)
  ie[[l]]<-rbind(est=as.matrix(temp2$ie)[l,],mean=apply(temp1$ie[[l]],2,mean),sd=apply(temp1$ie[[l]],2,sd),
                 upbd=apply(temp1$ie[[l]],2,mean)+b2*apply(temp1$ie[[l]],2,sd),
                 lwbd=apply(temp1$ie[[l]],2,mean)+b1*apply(temp1$ie[[l]],2,sd),
                 upbd_q=apply(temp1$ie[[l]],2,quantile,a2), lwbd_q=apply(temp1$ie[[l]],2,quantile,a1))
 names(ie)<-names(temp1$ie)
 temp1.result<-list(indirect.effect=ie,
                    total.effect=rbind(est=temp2$te,mean=apply(as.matrix(temp1$te),2,mean),sd=apply(as.matrix(temp1$te),2,sd),
                                   upbd=apply(as.matrix(temp1$te),2,mean)+b2*apply(as.matrix(temp1$te),2,sd),
                                   lwbd=apply(as.matrix(temp1$te),2,mean)+b1*apply(as.matrix(temp1$te),2,sd),
                                   upbd=apply(as.matrix(temp1$te),2,quantile,a2,na.rm=T),
                                   lwbd=apply(as.matrix(temp1$te),2,quantile,a1,na.rm=T)),
                    direct.effect=rbind(est=temp2$de,mean=apply(as.matrix(temp1$de),2,mean),sd=apply(as.matrix(temp1$de),2,sd),
                                   upbd=apply(as.matrix(temp1$de),2,mean)+b2*apply(as.matrix(temp1$de),2,sd),
                                   lwbd=apply(as.matrix(temp1$de),2,mean)+b1*apply(as.matrix(temp1$de),2,sd),
                                   upbd=apply(as.matrix(temp1$de),2,quantile,a2,na.rm=T),
                                   lwbd=apply(as.matrix(temp1$de),2,quantile,a1,na.rm=T)))
 ie<-NULL
 for (l in 1:nx)
   ie[[l]]<-rbind(est=as.matrix(temp4$ie)[l,],mean=apply(temp3$ie[[l]],2,mean),sd=apply(temp3$ie[[l]],2,sd),
                  upbd=apply(temp3$ie[[l]],2,mean)+b2*apply(temp3$ie[[l]],2,sd),
                  lwbd=apply(temp3$ie[[l]],2,mean)+b1*apply(temp3$ie[[l]],2,sd),
                  upbd_q=apply(temp3$ie[[l]],2,quantile,a2,na.rm=T), 
                  lwbd_q=apply(temp3$ie[[l]],2,quantile,a1,na.rm=T))
 names(ie)<-names(temp3$ie)
 temp2.result<-list(indirect.effect=ie,
                    direct.effect=rbind(est=temp4$de,mean=apply(as.matrix(temp3$de),2,mean),sd=apply(as.matrix(temp3$de),2,sd),
                                        upbd=apply(as.matrix(temp3$de),2,mean)+b2*apply(as.matrix(temp3$de),2,sd),
                                        lwbd=apply(as.matrix(temp3$de),2,mean)+b1*apply(as.matrix(temp3$de),2,sd),
                                        upbd=apply(as.matrix(temp3$de),2,quantile,a2,na.rm=T),
                                        lwbd=apply(as.matrix(temp3$de),2,quantile,a1,na.rm=T)))
 result<-list(results=temp1.result,re=temp2.result,alpha=alpha,plot=plot,obj=x,RE=RE,quant=quant,nx=nx,nie=nie,ny=ny)
 class(result)<-"summary.mma"
 result
 }

print.summary.mma<-function(x,...,digit=3)
{cat("MMA Analysis: Estimated Mediation Effects Using ")
 if (x$obj$model$MART)
  cat ("MART\n")
 else cat("GLM\n")
 pred.names<-names(x$result$indirect.effect)  
 for (l in 1:x$nx)
   {cat ("For Predictor",pred.names[l],"\n")
    temp.res<-list(total.effect=x$result$total.effect[,(x$ny*(l-1)+1):(x$ny*(l-1)+x$ny)],
                   direct.effect=x$result$direct.effect[,(x$ny*(l-1)+1):(x$ny*(l-1)+x$ny)],
                   indirect.effect=x$result$indirect.effect[[l]])
    print(lapply(temp.res,round,digit))
   }
 if(x$RE)
 {cat("The relative effects:\n")
   for (l in 1:x$nx)
   {cat ("For Predictor #",pred.names[l],"\n")
     temp.res<-list(direct.effect=x$re$direct.effect[,(x$ny*(l-1)+1):(x$ny*(l-1)+x$ny)],
                    indirect.effect=x$re$indirect.effect[[l]])
     print(lapply(temp.res,round,digit))
   }
 }
if(x$plot)
 for (l in 1:x$nx)
  for (m in 1:x$ny)
  {temp.t<-m%%x$ny
   temp.z<-(1:ncol(x$re$indirect.effect[[l]]))%%x$ny==temp.t
   temp.z[1:x$ny]<-F
   re<-c(x$re$indirect.effect[[l]][2,temp.z],x$re$dir[2,x$ny*(l-1)+m])
   if(x$quant)
    {upper<-c(x$re$indirect.effect[[l]][6,temp.z],x$re$dir[6,x$ny*(l-1)+m])
     lower<-c(x$re$indirect.effect[[l]][7,temp.z],x$re$dir[7,x$ny*(l-1)+m])}
   else
    {upper<-c(x$re$indirect.effect[[l]][4,temp.z],x$re$dir[4,x$ny*(l-1)+m])
     lower<-c(x$re$indirect.effect[[l]][5,temp.z],x$re$dir[5,x$ny*(l-1)+m])}
   d<-order(re)
   name1<-c(colnames(x$re$indirect.effect[[l]])[temp.z],"de")
   par(mfrow=c(1,1),mar=c(1,6,1,1),oma=c(3,2,2,4))
   bp <- barplot2(re[d], horiz = TRUE, main=paste("Relative Effects on y",m," on Predictor ",pred.names[l],sep=""), 
                names.arg=name1[d],plot.ci = TRUE, ci.u = upper[d], ci.l = lower[d],
                cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower)),
                col = rainbow(length(d), start = 3/6, end = 4/6))
}
}


plot.mma<-function(x,...,vari,xlim=range(x$data$x[,vari],na.rm=T),alpha=0.95,quantile=F)
{marg.den<-function(x,y,w=NULL) #added w
{if(!is.null(w))
  w<-w[!is.na(x) & !is.na(y)]
y<-y[!is.na(x)]
x<-x[!is.na(x)]
x<-x[!is.na(y)]
y<-y[!is.na(y)]
z1<-unique(x)
z2<-rep(0,length(z1))
if(is.null(w))   #
  for (i in 1:length(z1))
    z2[i]<-mean(y[x==z1[i]],na.rm=T)  
else          #
  for (i in 1:length(z1))      #
    z2[i]<-weighted.mean(y[x==z1[i]],w[x==z1[i]],na.rm=T)  #added ,w[x==z1[i]]
z3<-order(z1)
cbind(z1[z3],z2[z3])
}

weighted.hist<-function (x, w, breaks = "Sturges", col = NULL, plot = TRUE, 
                         freq = TRUE, ylim = NA, ylab = NULL, xaxis = TRUE, ...) 
{
  if (missing(x)) 
    stop("Usage: weighted.hist(x,...) vector of values x required")
  if (missing(w)) 
    w <- rep(1, length(x))
  breaks <- get.breaks(x, breaks)
  width <- diff(breaks)
  diffx <- diff(range(x))
  equidist <- sum(width - width[1]) < diffx/1000
  nbreaks <- length(breaks) - 1
  lastbreak <- breaks[nbreaks + 1]
  breaks[nbreaks + 1] <- breaks[nbreaks + 1] + diffx/1000
  if (diff(range(breaks)) < diffx) 
    warning("Not all values will be included in the histogram")
  counts <- rep(0, nbreaks)
  for (bin in 1:nbreaks) counts[bin] <- sum(w[x >= breaks[bin] & 
                                                x < breaks[bin + 1]])
  density <- counts/sum(counts)
  if (freq) {
    if (is.null(ylab)) 
      ylab <- "Frequency"
    heights <- counts
    if (!equidist) 
      warning("Areas will not relate to frequencies")
  }
  else {
    if (!equidist) {
      heights <- density * mean(width)/width
      heights <- heights/sum(heights)
    }
    else heights <- density
    if (is.null(ylab)) 
      ylab <- "Density"
  }
  if (plot) {
    if (is.null(col)) 
      col <- par("bg")
    if (is.na(ylim)) 
      ylim <- c(0, 1.1 * max(heights, na.rm = TRUE))
    mids <- barplot(heights, width = width, col = col, space = 0, 
                    ylim = ylim, ylab = ylab, ...)
    tickpos <- c(mids - width/2, mids[length(mids)] + width[length(width)]/2)
    if (xaxis) 
      axis(1, at = tickpos, labels = signif(c(breaks[1:nbreaks], 
                                              lastbreak), 3))
  }
  else mids <- breaks[-length(breaks)] + width/2
  invisible(list(breaks = breaks, counts = counts, density = density, 
                 mids = mids, xname = deparse(substitute(x)), equidist = equidist))
}

get.breaks<-function (x, breaks) 
{
  if (is.character(breaks)) 
    nbreaks <- do.call(paste("nclass", breaks, sep = ".", 
                             collapse = ""), list(x))
  if (is.numeric(breaks)) {
    if (length(breaks) == 1) {
      nbreaks <- breaks
    }
    else return(breaks)
  }
  breakinc <- diff(range(x))/nbreaks
  breaks <- c(min(x), rep(breakinc, nbreaks))
  breaks <- cumsum(breaks)
  return(breaks)
}

overlapHist <- function(a, b,breaks=NULL, xlim=NULL, xname=NULL, w=NULL)
{if(ncol(b)>1)
  {d<-rep(0,length(b))
   for (l in 1:ncol(b))
     d[b[,l]==1]<-l
   b<-d}
a1<-a
b1<-b
a<-a[!is.na(a1) & !is.na(b1)]
b<-b[!is.na(a1) & !is.na(b1)]
if(!is.null(w))                     #
  w<-w[!is.na(a1) & !is.na(b1)]    #
j<-sort(unique(b))
ahist<-hist(a[b==j[1]],plot=F)
if(!is.null(w))                     #
  ahist<-weighted.hist(a[b==j[1]], w[b==j[1]], plot=F)    #
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
if(is.null(w))                     #
  for (i in j)
    hist(a[b==i],ylab="Density",xlab="",breaks=breaks, 
         xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="="))
else           #
  for (i in j) #
    weighted.hist(a[b==i],w[b==i],ylab="Density",xlab="",breaks=breaks, #
                  xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="=")) #
}

weighted.prop.table<-function(x,w)  #the whole function is added for weighted proportions
{sumw<-sum(w)
temp<-sort(unique(x))
table<-c(0,length(temp))
names(table)<-temp
j<-1
for(temp1 in temp)
{table[j]<-sum(w[x==temp1])/sumw
j<-j+1}
table
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
{plot(df$x, df$F, ylim = range(c(df$L,df$U),na.rm=T), type = "l",xlab=xlab,ylab=ylab)
  polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
  lines(df$x, df$F, lwd = 2)
  lines(df$x, df$U, col="red",lty=2)
  lines(df$x, df$L, col="red",lty=2)}


nx<-ncol(x$data$dirx)
ny<-ncol(x$data$y)
op <- par(no.readonly = TRUE) # the whole list of settable par's.
data=x$data
mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
if (x$model[1]==T) 
 for (m in 1:ny) {
  full.model=x$model$model[[m]]
  best.iter=x$model$best.iter[m]
  if(data$binpred)
   {d<-rep(0,nrow(data$dirx))
    for(l in 1:nx)
      d[data$dirx[,l]==1]<-l
    if(!is.factor(data$x[,vari]))
     {par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      if(full.model$distribution=="gaussian")
        suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
      else if(full.model$distribution=="coxph")
        suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
      else
        suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
      overlapHist(a=data$x[,vari],b=as.matrix(d),xlim=xlim,xname="Predictor",w=data$w) # added w
     }
   else{par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    temp1<-NULL
    if (is.null(data$w)) #
    {temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])))
     for (j in 1:nx)
      temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1,vari])))
     barplot(prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])),ylim=c(0,max(temp1,na.rm=T)),
             ylab="Prop",sub=paste("Predictor at the Reference Level: pred=",0,sep=""))     
     for (j in 1:nx)
      barplot(prop.table(table(data$x[data$dirx[,j]==1,vari])),ylim=c(0,max(temp1,na.rm=T)),
             ylab="Prop",sub=paste(colnames(data$dirx)[j], ", pred=",j,sep=""))}
    else #
    {temp1<-c(temp1,weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari],data$w))#
     for (j in 1:nx) #
      temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w))#
     barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari]),ylim=c(0,max(temp1,na.rm=T)),#
             ylab="Prop",sub=paste("Predictor at the Reference Level, pred=", j,sep=""))
     for (j in 1:nx)#
      barplot(weighted.prop.table(data$x[data$dirx[,j]==1,vari]),ylim=c(0,max(temp1,na.rm=T)),#
              ylab="Prop",sub=paste(colnames(data$dirx)[j], ", pred=", j,sep=""))} #
  }
 }
else
{par(mfrow=c(3,nx),mar=c(5,5,1,1),oma=c(3,2,5,4))
  for (l in 1:nx)
   {temp.ie.detail<-as.matrix(x$boot.detail$ie1[[l]][,grep(mname,colnames(x$boot.detail$ie1[[l]]))])  #
    ie1<-boot.ci(x$boot.detail$pred.new[,l],matrix(temp.ie.detail[,m],nrow=nrow(x$boot.detail$pred.new)),alpha,quantile)
    plot_ci(ie1,xlab=colnames(data$dirx)[l],ylab=paste("IE on",colnames(data$y)[m]))}
  
  if(!is.factor(data$x[,vari]))
  {if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    if(nx>1)
      for (i in 1:(nx-1))
        plot(1, type="n", axes=F, xlab="", ylab="")
    for(l in 1:nx){
    axis(1,at=data$x[,vari],labels=F)
    a<-marg.den(data$dirx[,l],data$x[,vari],data$w) #added data$w
    scatter.smooth(a[,1],a[,2],family="gaussian",xlab=colnames(data$dirx)[l],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
  }
  else
  {if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    if(nx>1)
      for (i in 1:(nx-1))
        plot(1, type="n", axes=F, xlab="", ylab="")
    for(l in 1:nx){
      plot(data$x[,vari],data$dirx[,l],ylab=colnames(data$dirx)[l],xlab="")}}
}
}
else
  for (m in 1:ny) 
    {full.model=x$model$model[[m]]
     coef<-full.model$coefficients[names(full.model$coefficients)==vari] #plot the straight line instead of the loess line
     if(is.null(full.model$na.action))
       data1<-data$x[,vari]
     else
       data1<-data$x[-full.model$na.action,vari]
     if(data$binpred)
       {d<-rep(0,nrow(data$dirx))
        for(l in 1:nx)
         d[data$dirx[,l]==1]<-l
        if(!is.factor(data$x[,vari]))
        {par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
         if(!x$model$Survival[m])
            b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values),data$w) #added data$w
         else
            b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit,data$w)  #added data$w
         plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
         abline(a=mean(b[,2])-coef*mean(b[,1]),b=coef)
         axis(1,at=data1,labels=F)
         overlapHist(a=data$x[,vari],b=as.matrix(d),xlim=xlim,xname="Predictor",data$w)  #added data$w
        }
       else{par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
        if (!x$model$Survival[m])
          plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
        else
          plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
        temp1<-NULL
        if(is.null(data$w)){ #
           temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])))
           for (j in 1:ncol(data$dirx))
             temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1,vari])))
           barplot(prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])),ylim=c(0,max(temp1,na.rm=T)),
                   ylab="Prop",sub="Predictor at the reference level")
           for (j in 1:ncol(data$dirx))
             barplot(prop.table(table(data$x[data$dirx[,j]==1,vari])),ylim=c(0,max(temp1,na.rm=T)),
                     ylab="Prop",sub=colnames(data$dirx)[j])
       }#
       else#
          {temp1<-c(temp1,weighted.prop.table(table(data$x[apply(data$dirx,1,sum)==0,vari],data$w)))
           for (j in 1:ncol(data$dirx))#
             temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w))#
           barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari],data$w),ylim=c(0,max(temp1)),#
                   ylab="Prop",sub="Predictor at the reference level") #
           for (j in 1:ncol(data$dirx))#
             barplot(weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w),ylim=c(0,max(temp1)),#
                     ylab="Prop",sub=colnames(data$dirx)[j])} #
          }
      }
   else
    {par(mfrow=c(3,nx),mar=c(5,5,1,1),oma=c(3,2,5,4))
      for (l in 1:nx) {
         temp.ie.detail<-as.matrix(x$boot.detail$ie1[[l]][,grep(mname,colnames(x$boot.detail$ie1[[l]]))])  #
         ie1<-boot.ci(x$boot.detail$pred.new[,l],matrix(temp.ie.detail[,m],nrow=nrow(x$boot.detail$pred.new)),alpha,quantile)
         plot_ci(ie1,xlab=colnames(data$dirx)[l])}
      if(!is.factor(data$x[,vari]))
      {if(!x$model$Survival[m])
         b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values),data$w) #added data$w
       else
         b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit,data$w) #added data$w
       plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
       abline(a=mean(b[,2],na.rm=T)-coef*mean(b[,1],na.rm=T),b=coef)
       axis(1,at=data1,labels=F)
       if(nx>1)
         for (i in 1:(nx-1))
           plot(1, type="n", axes=F, xlab="", ylab="")
       for(l in 1:nx){
         a<-marg.den(data$dirx[,l],data$x[,vari],data$w)   #added data$w
         scatter.smooth(a[,1],a[,2],family="gaussian", xlab=colnames(data$dirx)[l],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
      }
    else
     {if (!x$model$Survival[m])
        plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
      else  
        plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data$x[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
      if(nx>1)
        for (i in 1:(nx-1))
          plot(1, type="n", axes=F, xlab="", ylab="")
      for(l in 1:nx){
         plot(data$x[,vari],data$dirx[,l],ylab=colnames(data$dirx)[l],xlab="")}}
}
}
par(op)
}



#plot on the med object
plot.med<-function(x,data,...,vari,xlim=range(data$x[,vari],na.rm=T))#data is the result from data.org
{marg.den<-function(x,y,w=NULL) #added w
{if(!is.null(w))
  w<-w[!is.na(x) & !is.na(y)]
y<-y[!is.na(x)]
x<-x[!is.na(x)]
x<-x[!is.na(y)]
y<-y[!is.na(y)]
z1<-unique(x)
z2<-rep(0,length(z1))
if(is.null(w))   #
  for (i in 1:length(z1))
    z2[i]<-mean(y[x==z1[i]],na.rm=T)  
else          #
  for (i in 1:length(z1))      #
    z2[i]<-weighted.mean(y[x==z1[i]],w[x==z1[i]],na.rm=T)  #added ,w[x==z1[i]]
z3<-order(z1)
cbind(z1[z3],z2[z3])
}

weighted.hist<-function (x, w, breaks = "Sturges", col = NULL, plot = TRUE, 
                         freq = TRUE, ylim = NA, ylab = NULL, xaxis = TRUE, ...) 
{
  if (missing(x)) 
    stop("Usage: weighted.hist(x,...) vector of values x required")
  if (missing(w)) 
    w <- rep(1, length(x))
  breaks <- get.breaks(x, breaks)
  width <- diff(breaks)
  diffx <- diff(range(x))
  equidist <- sum(width - width[1]) < diffx/1000
  nbreaks <- length(breaks) - 1
  lastbreak <- breaks[nbreaks + 1]
  breaks[nbreaks + 1] <- breaks[nbreaks + 1] + diffx/1000
  if (diff(range(breaks)) < diffx) 
    warning("Not all values will be included in the histogram")
  counts <- rep(0, nbreaks)
  for (bin in 1:nbreaks) counts[bin] <- sum(w[x >= breaks[bin] & 
                                                x < breaks[bin + 1]])
  density <- counts/sum(counts)
  if (freq) {
    if (is.null(ylab)) 
      ylab <- "Frequency"
    heights <- counts
    if (!equidist) 
      warning("Areas will not relate to frequencies")
  }
  else {
    if (!equidist) {
      heights <- density * mean(width)/width
      heights <- heights/sum(heights)
    }
    else heights <- density
    if (is.null(ylab)) 
      ylab <- "Density"
  }
  if (plot) {
    if (is.null(col)) 
      col <- par("bg")
    if (is.na(ylim)) 
      ylim <- c(0, 1.1 * max(heights, na.rm = TRUE))
    mids <- barplot(heights, width = width, col = col, space = 0, 
                    ylim = ylim, ylab = ylab, ...)
    tickpos <- c(mids - width/2, mids[length(mids)] + width[length(width)]/2)
    if (xaxis) 
      axis(1, at = tickpos, labels = signif(c(breaks[1:nbreaks], 
                                              lastbreak), 3))
  }
  else mids <- breaks[-length(breaks)] + width/2
  invisible(list(breaks = breaks, counts = counts, density = density, 
                 mids = mids, xname = deparse(substitute(x)), equidist = equidist))
}

get.breaks<-function (x, breaks) 
{
  if (is.character(breaks)) 
    nbreaks <- do.call(paste("nclass", breaks, sep = ".", 
                             collapse = ""), list(x))
  if (is.numeric(breaks)) {
    if (length(breaks) == 1) {
      nbreaks <- breaks
    }
    else return(breaks)
  }
  breakinc <- diff(range(x))/nbreaks
  breaks <- c(min(x), rep(breakinc, nbreaks))
  breaks <- cumsum(breaks)
  return(breaks)
}

overlapHist <- function(a, b,breaks=NULL, xlim=NULL, xname=NULL, w=NULL)
{if(ncol(b)>1)
{d<-rep(0,length(b))
for (l in 1:ncol(b))
  d[b[,l]==1]<-l
b<-d}
  a1<-a
  b1<-b
  a<-a[!is.na(a1) & !is.na(b1)]
  b<-b[!is.na(a1) & !is.na(b1)]
  if(!is.null(w))                     #
    w<-w[!is.na(a1) & !is.na(b1)]    #
  j<-sort(unique(b))
  ahist<-hist(a[b==j[1]],plot=F)
  if(!is.null(w))                     #
    ahist<-weighted.hist(a[b==j[1]], w[b==j[1]], plot=F)    #
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
  if(is.null(w))                     #
    for (i in j)
      hist(a[b==i],ylab="Density",xlab="",breaks=breaks, 
           xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="="))
  else           #
    for (i in j) #
      weighted.hist(a[b==i],w[b==i],ylab="Density",xlab="",breaks=breaks, #
                    xlim=xlim, ylim=c(0,yl), freq=F,main=paste(xname,i,sep="=")) #
}


weighted.prop.table<-function(x,w)  #the whole function is added for weighted proportions
{sumw<-sum(w)
temp<-sort(unique(x))
table<-c(0,length(temp))
names(table)<-temp
j<-1
for(temp1 in temp)
{table[j]<-sum(w[x==temp1])/sumw
j<-j+1}
table
}
op <- par(no.readonly = TRUE) # the whole list of settable par's.
nx<-length(x$ie)
ny<-ncol(data.frame(x$data$y))
data<-x$data
mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
if (x$model[1]==T) 
 for (m in 1:ny) {
  full.model=x$model$model[[m]]
  best.iter=x$model$best.iter[m]
  if(data$binpred)
  {d<-rep(0,nrow(data$dirx))
   for(l in 1:nx)
    d[data$dirx[,l]==1]<-l
   if(!is.factor(data$x[,vari]))
    {par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
     if(full.model$distribution=="gaussian")
        suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
     else if(full.model$distribution=="coxph")
        suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
     else
        suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
     overlapHist(a=data$x[,vari],b=as.matrix(d),xlim=xlim,xname="Predictor",w=data$w) #added w
    }
  else{par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(full.model$distribution=="gaussian")
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    temp1<-NULL
    if (is.null(data$w)) #
    {temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])))
    for (j in 1:nx)
      temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1,vari])))
    barplot(prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])),ylim=c(0,max(temp1,na.rm=T)),
            ylab="Prop",sub=paste("Predictor at the Reference Level: pred=",0,sep=""))     
    for (j in 1:nx)
      barplot(prop.table(table(data$x[data$dirx[,j]==1,vari])),ylim=c(0,max(temp1,na.rm=T)),
              ylab="Prop",sub=paste(colnames(data$dirx)[j], ", pred=",j,sep=""))}
    else #
    {temp1<-c(temp1,weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari],data$w))#
    for (j in 1:nx) #
      temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w))#
    barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari]),ylim=c(0,max(temp1,na.rm=T)),#
            ylab="Prop",sub=paste("Predictor at the Reference Level, pred=", j,sep=""))
    for (j in 1:nx)#
      barplot(weighted.prop.table(data$x[data$dirx[,j]==1,vari]),ylim=c(0,max(temp1,na.rm=T)),#
              ylab="Prop",sub=paste(colnames(data$dirx)[j], ", pred=", j,sep=""))} #
  }
}
else
{par(mfrow=c(3,nx),mar=c(5,5,1,1),oma=c(3,2,5,4)) #test
  for (l in 1:nx){
  temp2<-data$dirx[,l]
  temp3<-as.matrix(x$ie[[l]][,grep(mname,colnames(x$ie[[l]]))])
  temp3<-temp3[,m]
  temp.order=order(temp2)
  plot(temp2[temp.order], temp3[temp.order],type="l",
       xlab=colnames(data$dirx)[l],ylab=paste(c("IE of", vari, "on", 
       colnames(data$y)[m]),sep=""))}
  
  if(!is.factor(data$x[,vari]))
  {if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari,xlim=xlim))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,xlim=xlim,type="response"))
    if(nx>1)
      for (i in 1:(nx-1))
        plot(1, type="n", axes=F, xlab="", ylab="")
    for(l in 1:nx){
      axis(1,at=data$x[,vari],labels=F)
      a<-marg.den(data$dirx[,l],data$x[,vari],data$w) #added data$w
      scatter.smooth(a[,1],a[,2],family="gaussian",xlab=colnames(data$dirx)[l],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
  }
  else
  {if(full.model$distribution=="gaussian")
    suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter))
    else if(full.model$distribution=="coxph")
      suppressWarnings(plot.gbm(full.model, i.var=vari))
    else
      suppressWarnings(plot.gbm(full.model, i.var=vari,best.iter,type="response"))
    if(nx>1)
      for (i in 1:(nx-1))
        plot(1, type="n", axes=F, xlab="", ylab="")
    for(l in 1:nx){
      plot(data$x[,vari],data$dirx[,l],ylab=colnames(data$dirx)[l],xlab="")}}
}
}
else
  for (m in 1:ny) 
  {full.model=x$model$model[[m]]
   coef<-full.model$coefficients[names(full.model$coefficients)==vari] #plot the straight line instead of the loess line
   if(is.null(full.model$na.action))
     data1<-data$x[,vari]
   else
     data1<-data$x[-full.model$na.action,vari]
   if(data$binpred)
  {d<-rep(0,nrow(data$dirx))
   for(l in 1:nx)
    d[data$dirx[,l]==1]<-l
   if(!is.factor(data$x[,vari]))
    {par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
     if(!x$model$Survival[m])
      b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values),data$w) # added w
     else
      b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit,data$w) #added w
     plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
     abline(a=mean(b[,2])-coef*mean(b[,1]),b=coef)
     axis(1,at=data1,labels=F)
     overlapHist(a=data$x[,vari],b=as.matrix(d),xlim=xlim,xname="Predictor",data$w) #added w
   }
  else{par(mfrow=c(2+nx,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if (!x$model$Survival[m])
      plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    else
      plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    temp1<-NULL
    if(is.null(data$w)){ #
      temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])))
      for (j in 1:ncol(data$dirx))
        temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1,vari])))
      barplot(prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])),ylim=c(0,max(temp1,na.rm=T)),
              ylab="Prop",sub="Predictor at the reference level")
      for (j in 1:ncol(data$dirx))
        barplot(prop.table(table(data$x[data$dirx[,j]==1,vari])),ylim=c(0,max(temp1,na.rm=T)),
                ylab="Prop",sub=colnames(data$dirx)[j])
    }#
    else#
    {temp1<-c(temp1,weighted.prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari],data$w)))
    for (j in 1:ncol(data$dirx))#
      temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w))#
    barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0,vari],data$w),ylim=c(0,max(temp1)),#
            ylab="Prop",sub="Predictor at the reference level") #
    for (j in 1:ncol(data$dirx))#
      barplot(weighted.prop.table(data$x[data$dirx[,j]==1,vari],data$w),ylim=c(0,max(temp1)),#
              ylab="Prop",sub=colnames(data$dirx)[j])} #
  }
}
else
{par(mfrow=c(3,nx),mar=c(5,5,1,1),oma=c(3,2,5,4))
  for (l in 1:nx) {
    temp2<-data$dirx[,l]
    temp3<-as.matrix(x$ie[[l]][,grep(mname,colnames(x$ie[[l]]))])
    temp3<-temp3[,m]
    temp.order=order(temp2)
    plot(temp2[temp.order], temp3[temp.order],type="l",
         xlab=names(data$x)[data$dirx],ylab=paste(c("IE of", vari, "on", colnames(data$y)[m]),sep=""))}
  if(!is.factor(data$x[,vari]))
  {if(!x$model$Survival[m])
    b<-marg.den(data1,full.model$family$linkfun(full.model$fitted.values),data$w) #added data$w
   else
    b<-marg.den(data1,predict(full.model,se.fit=T,type=x$model$type)$fit,data$w) #added data$w
   plot(b,xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
   abline(a=mean(b[,2],na.rm=T)-coef*mean(b[,1],na.rm=T),b=coef)
   axis(1,at=data1,labels=F)
   if(nx>1)
     for (i in 1:(nx-1))
       plot(1, type="n", axes=F, xlab="", ylab="")
   for(l in 1:nx){
     a<-marg.den(data$x[,data$dirx],data$x[,vari],data$w) #added data$w
     scatter.smooth(a[,1],a[,2],family="gaussian", xlab=colnames(data$dirx)[l],ylim=xlim,ylab=paste("Mean",mname,sep="."))}
  }  
  else
  {if (!x$model$Survival[m])
    plot(full.model$fitted.values~data1,ylab=paste("f(",mname,")",sep=""),xlab=mname)
    else  
      plot(predict(full.model,se.fit=T,type=x$model$type)$fit~data$x[-full.model$na.action,vari],ylab=paste("f(",mname,")",sep=""),xlab=mname)
    if(nx>1)
      for (i in 1:(nx-1))
        plot(1, type="n", axes=F, xlab="", ylab="")
    for(l in 1:nx){
      plot(data$x[,vari],data$dirx[,l],ylab=colnames(data$dirx)[l],xlab="")}}
  }
}
par(op)
}


################################################################################################
# the functions to make inferences on mediation effects using parallel calculation
boot.med.par<-function(data,x=data$x, y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,catm=data$catm,
                       jointm=data$jointm,margin=1,n=20,seed=sample(1:1000,1),nonlinear=F,df=1,nu=0.001,
                       D=3,distn=NULL,family1=data$family1,n2=50,w=rep(1,nrow(x)),
                       refy=NULL,x.new=x,pred.new=dirx,binpred=data$binpred,type=NULL,w.new=NULL,ncore=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),
                         surv=(data$y_type==4),type,ncore=NULL)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, allm = c(contm, catm), 
                     n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn=NULL,family1=data$family1, #
                     biny=rep(F,ncol(y)),refy=rep(0,ncol(y)),surv=rep(F,ncol(y)),type=NULL, w=NULL) #
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    pred_names<-colnames(dirx)  #
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {te<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if (surv[m])
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      te
    }
    
    med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    new1<-nom1
    new1[,med]<-marg.m[1:nrow(nom1)]
    new0<-nom0
    new0[,med]<-marg.m[(nrow(nom1)+1):n3]
    dir.nom<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
    else if(surv[m])
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand)  
    {if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-as.factor(as.character(c(nom1[,med],nom0[,med])[temp.rand]))
    else
      marg.m<-c(nom1[,med],nom0[,med])[temp.rand]
    }        
      else                                         #added for the new program
        marg.m<-rbind(nom1[,med],nom0[,med])[temp.rand,]
      new1<-nom1
      new0<-nom0
      if(length(med)==1)                                       #added for the new program, when there is only one mediator
      {new1[,med]<-marg.m[1:nrow(new1)]                     #added for the new program 
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0))]}  #added for the new program
      else                                                     #added for the new program
      {new1[,med]<-marg.m[1:nrow(new1),]
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0)),]}
      dir.nom<-NULL
      for (m in 1:length(full.model))
        if(surv[m] & !is.null(best.iter1[m]))
          dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if(surv[m])
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      dir.nom
    }
    
    med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    temp.rand<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    marg.m1<-temp.rand[1:nrow(nom1)]
    marg.m2<-temp.rand[(nrow(nom1)+1):n3]
    dir.nom<-rep(0,length(full.model))
    for (m in 1:length(full.model))
      for (i in levels(x[,med]))
      {new1<-nom1
      new1[1:dim(new1)[1],med]<-i
      new0<-nom0
      new0[1:dim(new0)[1],med]<-i
      p<-mean(temp.rand==i,na.rm=T)
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
      else if(surv[m])
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T))
      else
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
      }
    dir.nom
    }
    
    #1.fit the model
    x2<-cbind(x,dirx)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    if (nonlinear)
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        full.model[[j]]<-suppressWarnings(gbm.fit(x1,y1, n.trees=200, interaction.depth=D, shrinkage=nu, w=w,
                                                  distribution=distn[j],train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))
        while(full.model[[j]]$n.trees-best.iter1[j]<30){
          full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
          best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}}
    else
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        if(surv[j])
          full.model[[j]]<-coxph(y1~., data=x1, weights=w)
        else
          full.model[[j]]<-glm(y1~., data=x1, family=family1[[j]], weights=w)
      }
    
    #2. prepare for the store of results
    set.seed(seed)
    te<-matrix(0,n,ncol(y)*ncol(dirx))
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    if(!is.null(jointm))
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
    }
    else
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)]),each=ncol(y)),sep=".")
    }
    denm<-rep(list(denm),ncol(dirx))
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x0.temp<-apply(dirx==1,1,sum)==0  #indicator of the reference group
      x0<-x2[x0.temp,]
      if(is.null(w))
      {w1<-NULL
      w0<-NULL}
      else
        w0<-w[x0.temp]
      for (l in 1:ncol(dirx))  #l indicate the lth predictor
      {x1.2<-x2[dirx[,l]==1,]
      if(!is.null(w))
        w1<-w[dirx[,l]==1]
      #n3<-dim(x)[1] use the original size
      new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
      new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
      te[k,((l-1)*ncol(y)+1):(l*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
      temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
      denm[[l]][k,1:ncol(y)]<-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand) #add temp.rand
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
        denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand)
        j<-j+1}
      #3.5 get the indirect effects
      ie[[l]][k,]<-te[k,((l-1)*ncol(y)+1):(l*ncol(y))]-denm[[l]][k,]
      if(!is.null(jointm))
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",names(x)[c(contm,catm)])
      }
    }
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),data=data)
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  if(is.null(jointm))
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))))
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie1)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  temp<-med.binx(data=NULL,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type,w=w)
  te[1,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie1[l,]<-apply(temp$ie[[l]],2,mean)}  #first row is the estimated value
  de[1,]<-apply(temp.1,2,mean,na.rm=T)
  model<-temp$model
  best.iter<-temp$best.iter
  
  registerDoParallel(cores=ncore)
  r<-foreach(i=1:n2,.combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T,prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  pred1<-data.frame(dirx[boots,])
  temp<-med.binx(data=NULL,x=x1, y=y1, dirx=pred1, contm=contm, catm=catm,jointm=jointm,allm=allm,n=n,seed=seed+i,
                 nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,biny=biny,refy=refy,surv=surv,type=type,w=NULL)
  te1<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  ie1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie1<-c(ie1,apply(temp$ie[[l]],2,mean,na.rm=T))}  #first row is the estimated value
  de1<-apply(temp.1,2,mean,na.rm=T)
  c(te1,de1,ie1)  
  }
  
  nxy=ny*nx
  nie=ncol(ie1)
  te[2:(n2+1),]<-r[,1:nxy]
  de[2:(n2+1),]<-r[,(nxy+1):(2*nxy)]
  temp.col<-2*nxy
  for (l in 1:nx)
  {ie[[l]][1:n2,]<-r[,(temp.col+1):(temp.col+nie)]
  temp.col<-temp.col+nie}
  
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         w=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type,w.new=NULL,ncore=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=data$family1,
                      biny=(data$y_type==2),refy=rep(NA,ncol(y)),x.new=x,pred.new=dirx, surv=(data$y_type==4),type=NULL,w=NULL, w.new=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    pred_names<-colnames(dirx)
    ynames<-colnames(y)
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
    
    anymissing<-function(vec) #return T if there is any missing in the vec
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row,w=NULL)
    {temp<-matrix(col,n.row)
    if(is.null(w))
      return(apply(temp,1,mean,na.rm=T))
    else
      return(apply(temp,1,weighted.mean,na.rm=T,w=w))}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df,w) #give the model and residual of m given x
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=binomial(link = "logit"),weights=w)
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=gaussian(link="identity"),weights=w)
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df)  #added nonlinear and df to sim.xm
    {mult.norm<-function(mu,vari,n) 
    {if (nrow(vari)!=ncol(vari)) 
      result<-c("Error: Variance matrix is not square")  
    else if (length(mu)!=nrow(vari)) 
      result<-c("Error: length mu is not right!")  
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
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    
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
    
    if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
    tempm[[1]]<-NULL}
    else  tempm<-NULL
    allm<-unique(c(contm,binm,unlist(tempm)))
    
    nonmissing<-apply(cbind(y,x[,listm$single],dirx),1,anymissing)
    x<-x[nonmissing,]
    y<-data.frame(y[nonmissing,])
    colnames(y)<-ynames
    pred<-data.frame(dirx[nonmissing,])
    colnames(pred)<-pred_names
    w<-w[nonmissing]
    nonmissing1<-apply(cbind(x.new[,listm$single],pred.new),1,anymissing)
    x.new<-x.new[nonmissing1,]
    w.new<-w.new[nonmissing1]
    pred.new<-data.frame(pred.new[nonmissing1,])
    colnames(pred.new)<-pred_names
    
    #1.fit the model
    x2<-cbind(x,pred)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    for(j in 1:ncol(y)){
      if(biny[j])                     #recode y if y is binary
        y[,j]<-ifelse(y[,j]==refy[j],0,1)
      
      if(nonlinear)
      {full.model[[j]]<-suppressWarnings(gbm.fit(x2,y[,j], n.trees=200, interaction.depth=D, shrinkage=nu,w=w,
                                                 distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
      best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))         
      while(full.model[[j]]$n.trees-best.iter1[j]<30){
        full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}
      }
      else
      {if(surv[j])
        full.model[[j]]<-coxph(y[,j]~., data=x2, weights=w)
      else
        full.model[[j]]<-glm(y[,j]~., data=x2, family=family1[[j]], weights=w)
      }
    }
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    #  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df,w)
    te1.0<-NULL
    denm1.0<-NULL
    
    n1<-dim(x)[1]
    
    #4. repeat to get the mediation effect
    for (l in 1:ncol(pred)) {
      denm1<-NULL
      te1<-NULL
      for (k in 1:n)
      {new0<-sim.xm(distmgivenx,x.new,pred.new,binm,contm,catm,nonlinear,df) #draw ms conditional on x.new
      temp.pred<-pred.new
      temp.pred[,l]<-temp.pred[,l]+margin
      new1<-sim.xm(distmgivenx,x.new,temp.pred,binm,contm,catm,nonlinear,df)  #draw from the conditional distribution of m given x
      new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
      new0<-cbind(new0,pred.new) 
      denm2<-NULL
      
      sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new)   #random sample from the original data
      
      temp.new1<-new1
      temp.new1[,allm]<-x.new[sample.temp[1:n.new],allm]
      temp.new0<-new0
      temp.new0[,allm]<-x.new[sample.temp[(n.new+1):(2*n.new)],allm]
      #4.0 get the direct effect
      for (m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],temp.new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)
      
      #4.1 get the te
      te0<-NULL
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
      te1<-cbind(te1,te0)
      
      #4.2 mediation effect from the single mediator
      if (!is.null(listm$single))
        for (i in 1:length(listm$single))
        {new1.nm<-new1
        new0.nm<-new0
        temp.m<-x.new[sample.temp,listm$single[i]]
        new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
        new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
        for(m in 1:ncol(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      
      #4.3.mediation effect from the joint mediator
      if (!is.null(listm$multi))
        for (i in 2:(listm$multi[[1]]+1))
        {new1.nm<-new1
        new0.nm<-new0
        new1.nm[,listm$multi[[i]]]<-x.new[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
        new0.nm[,listm$multi[[i]]]<-x.new[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
        for(m in 1:col(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      denm1<-rbind(denm1,denm2)
      }
      denm1.0[[l]]<-denm1 
      te1.0[[l]]<-te1
    } 
    
    #4.4 get the indirect effects
    denm<-NULL
    te<-NULL
    ie<-NULL
    for (l in 1:ncol(pred))
    {denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
    te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
    te<-cbind(te,te0)
    temp1<-ncol(denm[[l]])/ncol(te0)
    temp2<-NULL
    for(temp in 1:temp1)
      temp2<-cbind(temp2,te0)
    ie[[l]]<-temp2-denm[[l]]
    if(!is.null(listm$multi)) 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single]),each=ncol(y)),sep=".")
    if(!is.null(listm$multi))
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single]),each=ncol(y)),sep=".")
    }
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),pred.new=pred.new,w.new=w.new,data=data)
    class(a)<-"med"
    return(a)
  }
  
  anymissing<-function(vec)
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
  col_mean1<-function(col,n.row,w=NULL)
  {temp<-matrix(col,n.row)
  if(is.null(w))
    return(apply(temp,2,mean,na.rm=T))
  else
    return(apply(temp,2,weighted.mean,na.rm=T,w=w))}
  
  stock<-function(mat,rows)
  {if(nrow(mat)==rows)
    temp=mat
  else{
    temp<-mat[1:rows,]
    temp1=rows
    while(temp1<nrow(mat))   
    {temp<-cbind(temp,mat[(temp1+1):(temp1+rows),])
      temp1<-temp1+rows}
  }
  temp
  }
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
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
  
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2,ny*(1+length(listm$single)+mul))   #added in the new program
  ie1<-matrix(0,nx,ny*(1+length(listm$single)+mul))   #added in the new program
  if(!is.null(listm$multi))
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  
  temp<-med.contx(data=NULL,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,pred.new=pred.new, surv,type,w=w,w.new=w.new)
  
  temp.1<-NULL
  for (l in 1:nx)
    temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  if(is.null(w.new))
  {te[1,]<-apply(temp$te,2,mean,na.rm=T)
  de[1,]<-apply(temp.1,2,mean,na.rm=T) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
  }
  else
  {te[1,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
  de[1,]<-apply(temp$denm[,1],2,weighted.mean,na.rm=T,w=w.new) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new)  #first row is the estimated value
  }
  
  te1<-NULL                      #to store the mediation effects on predictor
  de1<-NULL
  ie2<-rep(list(NULL),nx)
  names(ie2)<-pred_names
  model<-temp$model
  
  registerDoParallel(cores=ncore)
  l<-1
  r<-foreach(l=1:n2, .combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T, prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  dirx1<-data.frame(dirx[boots,])
  temp<-med.contx(data=NULL,x=x1,y=y1,dirx=dirx1,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed+i,nonlinear=nonlinear,df=df,nu=nu,D=D,
                  distn=distn,family1=family1,biny=biny,refy=refy,x.new=x.new,pred.new=pred.new,surv=surv,type=type) #added to the new codel, change the seed to make different results
  temp.1<-NULL
  temp.2<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  temp.2<-cbind(temp.2,temp$ie[[l]])}
  te1<-temp$te  #n*nxy
  de1<-temp.1   #n*nxy
  ie1<-temp.2   #n*(nxy*nie)
  cbind(te1,de1,ie1)
  }
  nxy=nx*ny
  te1<-stock(as.matrix(r[,1:nxy]),nrow(x.new))
  de1<-stock(as.matrix(r[,(1+nxy):(2*nxy)]),nrow(x.new))
  nie<-ncol(ie[[1]])
  temp.count<-2*nxy
  for(l in 1:nx)
  {ie2[[l]]<-r[,(temp.count+1):(temp.count+nie)]
  temp.count=temp.count+nie}
  
  if(is.null(w.new))
  {te[2:(n2+1),]<-matrix(apply(te1,2,mean,na.rm=T),ncol=nxy,byrow=T)
  de[2:(n2+1),]<-matrix(apply(de1,2,mean,na.rm=T),ncol=nxy,byrow=T)
  for (l in 1:nx)
    ie[[l]]<-apply(ie2[[l]],2,col_mean1,nrow(x.new))  #first row is the estimated value
  }
  else
  {te[2:(n2+1),]<-matrix(apply(te1,2,weighted.mean,na.rm=T,w=w.new),ncol=nxy,byrow=T)
  de[2:(n2+1),]<-matrix(apply(de1,2,weighted.mean,na.rm=T,w=w.new),ncol=nxy,byrow=T)
  for (l in 1:nx)
    ie[[l]]<-apply(ie2[[l]],2,col_mean1,nrow(x.new),w=w.new)  #first row is the estimated value
  }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  missing.pred.new<-apply(pred.new,1,anymissing)
  pred.new<-data.frame(pred.new[!missing.pred.new,])
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(pred.new=pred.new,te1=te1,de1=de1,ie1=ie2),w.new=w.new)
  class(a)<-"mma"
  return(a)
}

if(is.null(data)){
  surv=rep(F,ncol(y))
  biny=rep(F,ncol(y))
  if(is.null(distn))
    distn<-rep(NA,ncol(y))
  for(j in 1:ncol(y)) {
    if(class(y[,j])=="Surv"){
      surv[j]=T
      if(is.na(distn[j]))
        distn[j]="coxph"
      if(is.null(type) & nonlinear)
        type="response"
      else if (is.null(type))
        type="risk"
    }
    else if(is.character(y[,j]) | is.factor(y[,j]) | nlevels(as.factor(y[,j]))==2)
    {biny[j]=T
    if(is.na(family1[[j]]))
      family1[[j]] = binomial("logit")
    if(is.na(distn[j]))
      distn[j]="bernoulli" 
    if(!is.na(refy[j]))
      y[,j]<-ifelse(y[,j]==refy[j],0,1)
    else
      y[,j]<-ifelse(as.factor(y[,j])==levels(as.factor(y[,j]))[1],0,1)
    }
    else { 
      if(is.na(family1[[j]]))
        family1[[j]] = gaussian(link = "identity")
      if(is.na(distn[j]))
        distn[j]="gaussian" 
    }
  }
}
else
{biny=data$y_type==2
surv=data$y_type==4
if(sum(surv)>0 & is.null(type) & nonlinear)
  type="response"
else if (sum(surv)>0 & is.null(type))
  type="risk"
if(is.null(distn))
  distn<-rep(NA,ncol(y))
distn[is.na(distn) & data$y_type==2]="bernoulli"
distn[is.na(distn) & data$y_type==4]="coxph"
distn[is.na(distn) & data$y_type==1]="gaussian"
}

if(binpred)
  a<-boot.med.binx(data=data,x=x, y=y,dirx=dirx,contm=contm,catm=catm,
                   jointm=jointm,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,
                   D=D,distn=distn,family1=family1,
                   w=w,biny=biny,refy=rep(0,ncol(y)),surv,type)
else
  a<-boot.med.contx(data=data,x=x,y=y,dirx=dirx,binm=binm,contm=contm,
                    catm=catm, jointm=jointm, margin = margin, n = n, seed = seed, 
                    nonlinear = nonlinear, df = df, nu = nu, D = D, distn = distn, 
                    family1 = family1, n2 = n2,w=w,
                    biny=biny,refy=rep(0,ncol(y)),x.new=x.new,pred.new=pred.new, surv,type,w.new)

return(a)
}


mma.par<-function(x,y,pred,mediator=NULL, contmed=NULL,binmed=NULL,binref=NULL,
                  catmed=NULL,catref=NULL,jointm=NULL,refy=rep(NA,ncol(data.frame(y))),
                  predref=NULL,alpha=0.1,alpha2=0.1, margin=1, n=20,seed=sample(1:1000,1),
                  nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=as.list(rep(NA,ncol(data.frame(y)))),
                  n2=50,w=rep(1,nrow(x)), testtype=1, x.new=NULL, pred.new=NULL, type=NULL,w.new=NULL, ncore=NULL)
{boot.med.binx<-function(data,x=data$x, y=data$y,dirx=data$dirx,contm=data$contm,catm=data$catm,
                         jointm=data$jointm,n=20,seed=sample(1:1000,1),n2=50,nonlinear=F,nu=0.001,
                         D=3,distn="bernoulli",family1=binomial("logit"),
                         w=rep(1,nrow(x)),biny=(data$y_type==2),refy=rep(NA,ncol(y)),
                         surv=(data$y_type==4),type,ncore=NULL)
  #n2 is the time of bootstrap
{
  med.binx<-function(data, x=data$x, y=data$y, dirx=data$dirx, contm = data$contm, 
                     catm = data$catm, jointm = data$jointm, allm = c(contm, catm), 
                     n=20,seed=sample(1:1000,1),nonlinear=F,nu=0.001,
                     D=3,distn=NULL,family1=data$family1, #
                     biny=rep(F,ncol(y)),refy=rep(0,ncol(y)),surv=rep(F,ncol(y)),type=NULL, w=NULL) #
  {if (is.null(allm))
    stop("Error: no potential mediator is specified")
    xnames<-colnames(x)
    pred_names<-colnames(dirx)  #
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(catm))
      catm<-unlist(sapply(catm,grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    allm=c(contm,catm)
    
    te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
    {te<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if (surv[m])
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      te
    }
    
    med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    new1<-nom1
    new1[,med]<-marg.m[1:nrow(nom1)]
    new0<-nom0
    new0[,med]<-marg.m[(nrow(nom1)+1):n3]
    dir.nom<-NULL
    for(m in 1:length(full.model))
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
    else if(surv[m])
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
    else
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
    dir.nom
    }
    
    med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand)  
    {if (length(med)==1)                       #added for the new program, when there is only one mediator
    {if(is.factor(nom1[,med]))              #added to control for one factor mediator
      marg.m<-as.factor(as.character(c(nom1[,med],nom0[,med])[temp.rand]))
    else
      marg.m<-c(nom1[,med],nom0[,med])[temp.rand]
    }        
      else                                         #added for the new program
        marg.m<-rbind(nom1[,med],nom0[,med])[temp.rand,]
      new1<-nom1
      new0<-nom0
      if(length(med)==1)                                       #added for the new program, when there is only one mediator
      {new1[,med]<-marg.m[1:nrow(new1)]                     #added for the new program 
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0))]}  #added for the new program
      else                                                     #added for the new program
      {new1[,med]<-marg.m[1:nrow(new1),]
      new0[,med]<-marg.m[(nrow(new1)+1):(nrow(new1)+nrow(new0)),]}
      dir.nom<-NULL
      for (m in 1:length(full.model))
        if(surv[m] & !is.null(best.iter1[m]))
          dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
      else if(surv[m])
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)
      else
        dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
      dir.nom
    }
    
    med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type)  
    {n3<-nrow(nom1)+nrow(nom0)
    temp.rand<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
    marg.m1<-temp.rand[1:nrow(nom1)]
    marg.m2<-temp.rand[(nrow(nom1)+1):n3]
    dir.nom<-rep(0,length(full.model))
    for (m in 1:length(full.model))
      for (i in levels(x[,med]))
      {new1<-nom1
      new1[1:dim(new1)[1],med]<-i
      new0<-nom0
      new0[1:dim(new0)[1],med]<-i
      p<-mean(temp.rand==i,na.rm=T)
      if(surv[m] & !is.null(best.iter1[m]))
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
      else if(surv[m])
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit,na.rm=T))
      else
        dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
      }
    dir.nom
    }
    
    #1.fit the model
    x2<-cbind(x,dirx)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    if (nonlinear)
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        full.model[[j]]<-suppressWarnings(gbm.fit(x1,y1, n.trees=200, interaction.depth=D, shrinkage=nu, w=w,
                                                  distribution=distn[j],train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))
        while(full.model[[j]]$n.trees-best.iter1[j]<30){
          full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
          best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}}
    else
      for (j in 1:ncol(y)){
        if(biny[j])                     #recode y if y is binary
          y[,j]<-ifelse(y[,j]==refy[j],0,1)
        x1<-x2[!is.na(y[,j]),]             #delete nas in y for mart
        y1<-y[!is.na(y[,j]),j]
        if(surv[j])
          full.model[[j]]<-coxph(y1~., data=x1, weights=w)
        else
          full.model[[j]]<-glm(y1~., data=x1, family=family1[[j]], weights=w)
      }
    
    #2. prepare for the store of results
    set.seed(seed)
    te<-matrix(0,n,ncol(y)*ncol(dirx))
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    if(!is.null(jointm))
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
    }
    else
    {denm<-matrix(0,n,ncol(y)*(1+length(c(contm,catm))))
    dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[c(contm,catm)]),each=ncol(y)),sep=".")
    }
    denm<-rep(list(denm),ncol(dirx))
    ie<-denm
    #3. repeat to get the mediation effect
    for (k in 1:n)
    {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
      x0.temp<-apply(dirx==1,1,sum)==0  #indicator of the reference group
      x0<-x2[x0.temp,]
      if(is.null(w))
      {w1<-NULL
      w0<-NULL}
      else
        w0<-w[x0.temp]
      for (l in 1:ncol(dirx))  #l indicate the lth predictor
      {x1.2<-x2[dirx[,l]==1,]
      if(!is.null(w))
        w1<-w[dirx[,l]==1]
      #n3<-dim(x)[1] use the original size
      new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
      new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
      te[k,((l-1)*ncol(y)+1):(l*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
      temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
      denm[[l]][k,1:ncol(y)]<-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand) #add temp.rand
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
        denm[[l]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand)
        j<-j+1}
      #3.5 get the indirect effects
      ie[[l]][k,]<-te[k,((l-1)*ncol(y)+1):(l*ncol(y))]-denm[[l]][k,]
      if(!is.null(jointm))
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie[[l]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",names(x)[c(contm,catm)])
      }
    }
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear, Survival=surv, type=type, model=full.model,best.iter=best.iter1),data=data)
    class(a)<-"med"
    return(a)
  }
  
  if (is.null(c(contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  if(is.character(contm))
    contm<-unlist(sapply(contm,grep,xnames))
  if(is.character(catm))
    catm<-unlist(sapply(catm,grep,xnames))
  if(!is.null(jointm))
    for (i in 2:length(jointm))
      if(is.character(jointm[[i]]))
        jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
  
  allm=c(contm,catm)
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  if(is.null(jointm))
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))))
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {ie<-matrix(0,n2,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  ie1<-matrix(0,nx,ny*(1+length(c(contm,catm))+jointm[[1]]))
  dimnames(ie1)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  temp<-med.binx(data=NULL,x,y,dirx,contm,catm,jointm,allm,n,seed,nonlinear,nu,D,distn,family1,biny,refy,surv,type,w=w)
  te[1,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie1[l,]<-apply(temp$ie[[l]],2,mean)}  #first row is the estimated value
  de[1,]<-apply(temp.1,2,mean,na.rm=T)
  model<-temp$model
  best.iter<-temp$best.iter
  
  registerDoParallel(cores=ncore)
  r<-foreach(i=1:n2,.combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T,prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  pred1<-data.frame(dirx[boots,])
  temp<-med.binx(data=NULL,x=x1, y=y1, dirx=pred1, contm=contm, catm=catm,jointm=jointm,allm=allm,n=n,seed=seed+i,
                 nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,biny=biny,refy=refy,surv=surv,type=type,w=NULL)
  te1<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-NULL
  ie1<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  ie1<-c(ie1,apply(temp$ie[[l]],2,mean,na.rm=T))}  #first row is the estimated value
  de1<-apply(temp.1,2,mean,na.rm=T)
  c(te1,de1,ie1)  
  }
  
  nxy=ny*nx
  nie=ncol(ie1)
  te[2:(n2+1),]<-r[,1:nxy]
  de[2:(n2+1),]<-r[,(nxy+1):(2*nxy)]
  temp.col<-2*nxy
  for (l in 1:nx)
  {ie[[l]][1:n2,]<-r[,(temp.col+1):(temp.col+nie)]
  temp.col<-temp.col+nie}
  
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model, 
          data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T))
  class(a)<-"mma"
  return(a)
}

boot.med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                         catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                         nonlinear=F,df=1,nu=0.001,D=3,distn="gaussian",
                         family1=gaussian(link="identity"),n2=50,
                         w=rep(1,nrow(x)),biny=F,refy=0,x.new=x,surv,type,w.new=NULL,ncore=NULL)
{
  med.contx<-function(data,x=data$x,y=data$y,dirx=data$dirx,binm=data$binm,contm=data$contm,
                      catm=data$catm, jointm=data$jointm, margin=1, n=20,seed=sample(1:1000,1),
                      nonlinear=F,df=1,nu=0.001,D=3,distn=NULL,family1=data$family1,
                      biny=(data$y_type==2),refy=rep(NA,ncol(y)),x.new=x,pred.new=dirx, surv=(data$y_type==4),type=NULL,w=NULL, w.new=NULL)
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
    
    xnames<-colnames(x)
    pred_names<-colnames(dirx)
    ynames<-colnames(y)
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
    
    anymissing<-function(vec) #return T if there is any missing in the vec
    {if(sum(is.na(vec))>0)
      return(F)
      else return(T)
    }
    
    col_mean<-function(col,n.row,w=NULL)
    {temp<-matrix(col,n.row)
    if(is.null(w))
      return(apply(temp,1,mean,na.rm=T))
    else
      return(apply(temp,1,weighted.mean,na.rm=T,w=w))}
    
    
    dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df,w) #give the model and residual of m given x
    {models<-NULL
    res<-NULL
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    j<-1
    if(!is.null(binm))
    {for(i in binm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=binomial(link = "logit"),weights=w)
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
    j<-j+1}
    }
    for (i in contm)
    { models[[j]]<-glm(x[,i]~.,data=data.frame(z),family=gaussian(link="identity"),weights=w)
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
    }
    
    
    sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df)  #added nonlinear and df to sim.xm
    {mult.norm<-function(mu,vari,n) 
    {if (nrow(vari)!=ncol(vari)) 
      result<-c("Error: Variance matrix is not square")  
    else if (length(mu)!=nrow(vari)) 
      result<-c("Error: length mu is not right!")  
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
    if(nonlinear)
    {z<-NULL
    for(i in 1:ncol(dirx))
      z<-cbind(z,ns(dirx[,i],df=df))}
    else
      z<-dirx
    
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
    
    if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
    tempm[[1]]<-NULL}
    else  tempm<-NULL
    allm<-unique(c(contm,binm,unlist(tempm)))
    
    nonmissing<-apply(cbind(y,x[,listm$single],dirx),1,anymissing)
    x<-x[nonmissing,]
    y<-data.frame(y[nonmissing,])
    colnames(y)<-ynames
    pred<-data.frame(dirx[nonmissing,])
    colnames(pred)<-pred_names
    w<-w[nonmissing]
    nonmissing1<-apply(cbind(x.new[,listm$single],pred.new),1,anymissing)
    x.new<-x.new[nonmissing1,]
    w.new<-w.new[nonmissing1]
    pred.new<-data.frame(pred.new[nonmissing1,])
    colnames(pred.new)<-pred_names
    
    #1.fit the model
    x2<-cbind(x,pred)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-NULL
    best.iter1<-NULL
    for(j in 1:ncol(y)){
      if(biny[j])                     #recode y if y is binary
        y[,j]<-ifelse(y[,j]==refy[j],0,1)
      
      if(nonlinear)
      {full.model[[j]]<-suppressWarnings(gbm.fit(x2,y[,j], n.trees=200, interaction.depth=D, shrinkage=nu,w=w,
                                                 distribution=distn,train.fraction=1.0, bag.fraction=0.5, verbose=FALSE))
      best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))         
      while(full.model[[j]]$n.trees-best.iter1[j]<30){
        full.model[[j]]<-gbm.more(full.model[[j]], 50)           # do another 50 iterations
        best.iter1[j]<-suppressWarnings(gbm.perf(full.model[[j]],plot.it=FALSE,method="OOB"))}
      }
      else
      {if(surv[j])
        full.model[[j]]<-coxph(y[,j]~., data=x2, weights=w)
      else
        full.model[[j]]<-glm(y[,j]~., data=x2, family=family1[[j]], weights=w)
      }
    }
    
    #2. prepare for the store of results
    set.seed(seed)
    n.new<-nrow(x.new)
    #  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
    
    #3. get the joint distribution of m given x
    distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df,w)
    te1.0<-NULL
    denm1.0<-NULL
    
    n1<-dim(x)[1]
    
    #4. repeat to get the mediation effect
    for (l in 1:ncol(pred)) {
      denm1<-NULL
      te1<-NULL
      for (k in 1:n)
      {new0<-sim.xm(distmgivenx,x.new,pred.new,binm,contm,catm,nonlinear,df) #draw ms conditional on x.new
      temp.pred<-pred.new
      temp.pred[,l]<-temp.pred[,l]+margin
      new1<-sim.xm(distmgivenx,x.new,temp.pred,binm,contm,catm,nonlinear,df)  #draw from the conditional distribution of m given x
      new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
      new0<-cbind(new0,pred.new) 
      denm2<-NULL
      
      sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new)   #random sample from the original data
      
      temp.new1<-new1
      temp.new1[,allm]<-x.new[sample.temp[1:n.new],allm]
      temp.new0<-new0
      temp.new0[,allm]<-x.new[sample.temp[(n.new+1):(2*n.new)],allm]
      #4.0 get the direct effect
      for (m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],temp.new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)
      
      #4.1 get the te
      te0<-NULL
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
      else if(surv[m])
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
      else
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
      te1<-cbind(te1,te0)
      
      #4.2 mediation effect from the single mediator
      if (!is.null(listm$single))
        for (i in 1:length(listm$single))
        {new1.nm<-new1
        new0.nm<-new0
        temp.m<-x.new[sample.temp,listm$single[i]]
        new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
        new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
        for(m in 1:ncol(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      
      #4.3.mediation effect from the joint mediator
      if (!is.null(listm$multi))
        for (i in 2:(listm$multi[[1]]+1))
        {new1.nm<-new1
        new0.nm<-new0
        new1.nm[,listm$multi[[i]]]<-x.new[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
        new0.nm[,listm$multi[[i]]]<-x.new[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
        for(m in 1:col(y))
          if(surv[m] & !is.null(best.iter1[m]))
            denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
        else if(surv[m])
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type,se.fit=TRUE)$fit-predict(full.model[[m]],new0.nm,best.iter1[m],type=type,se.fit=TRUE)$fit)/margin)
        else
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
        }
      denm1<-rbind(denm1,denm2)
      }
      denm1.0[[l]]<-denm1 
      te1.0[[l]]<-te1
    } 
    
    #4.4 get the indirect effects
    denm<-NULL
    te<-NULL
    ie<-NULL
    for (l in 1:ncol(pred))
    {denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
    te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
    te<-cbind(te,te0)
    temp1<-ncol(denm[[l]])/ncol(te0)
    temp2<-NULL
    for(temp in 1:temp1)
      temp2<-cbind(temp2,te0)
    ie[[l]]<-temp2-denm[[l]]
    if(!is.null(listm$multi)) 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",names(x)[listm$single]),each=ncol(y)),sep=".")
    if(!is.null(listm$multi))
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
    else 
      colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",names(x)[listm$single]),each=ncol(y)),sep=".")
    }
    colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
    names(denm)<-pred_names
    names(ie)<-pred_names
    a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, model=full.model,best.iter=best.iter1),pred.new=pred.new,w.new=w.new,data=data)
    class(a)<-"med"
    return(a)
  }
  
  anymissing<-function(vec)
  {if(sum(is.na(vec))>0)
    return(F)
    else return(T)
  }
  
  col_mean1<-function(col,n.row,w=NULL)
  {temp<-matrix(col,n.row)
  if(is.null(w))
    return(apply(temp,2,mean,na.rm=T))
  else
    return(apply(temp,2,weighted.mean,na.rm=T,w=w))}
  
  stock<-function(mat,rows)
  {if(nrow(mat)==rows)
    temp=mat
  else{
    temp<-mat[1:rows,]
    temp1=rows
    while(temp1<nrow(mat))   
    {temp<-cbind(temp,mat[(temp1+1):(temp1+rows),])
      temp1<-temp1+rows}
  }
  temp
  }
  
  if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
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
  
  ny=ncol(y)
  nx=ncol(dirx)
  te<-matrix(0,n2+1,ny*nx)
  de<-matrix(0,n2+1,ny*nx)
  mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
  ie<-matrix(0,n2,ny*(1+length(listm$single)+mul))   #added in the new program
  ie1<-matrix(0,nx,ny*(1+length(listm$single)+mul))   #added in the new program
  if(!is.null(listm$multi))
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single],name1),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  else 
  {dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",names(x)[listm$single]),each=ny),sep=".")
  rownames(ie1)<-pred_names}
  ie<-rep(list(ie),nx)
  names(ie)<-pred_names
  
  temp<-med.contx(data=NULL,x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed,nonlinear=nonlinear,df=df,nu=nu,D=D,distn=distn,family1=family1,biny=biny,
                  refy=refy,x.new=x.new,pred.new=pred.new, surv,type,w=w,w.new=w.new)
  
  temp.1<-NULL
  for (l in 1:nx)
    temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  if(is.null(w.new))
  {te[1,]<-apply(temp$te,2,mean,na.rm=T)
  de[1,]<-apply(temp.1,2,mean,na.rm=T) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
  }
  else
  {te[1,]<-apply(temp$te,2,weighted.mean,na.rm=T,w=w.new)
  de[1,]<-apply(temp$denm[,1],2,weighted.mean,na.rm=T,w=w.new) 
  for (l in 1:nx)
    ie1[l,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new)  #first row is the estimated value
  }
  
  te1<-NULL                      #to store the mediation effects on predictor
  de1<-NULL
  ie2<-rep(list(NULL),nx)
  names(ie2)<-pred_names
  model<-temp$model
  
  registerDoParallel(cores=ncore)
  l<-1
  r<-foreach(l=1:n2, .combine=rbind,.packages =c('gbm','survival','splines','doParallel'))  %dopar% 
  {boots<-sample(1:nrow(x),replace=T, prob=w)
  x1<-x[boots,]
  y1<-data.frame(y[boots,])
  dirx1<-data.frame(dirx[boots,])
  temp<-med.contx(data=NULL,x=x1,y=y1,dirx=dirx1,binm=binm,contm=contm,catm=catm,jointm=jointm, 
                  margin=margin,n=n,seed=seed+i,nonlinear=nonlinear,df=df,nu=nu,D=D,
                  distn=distn,family1=family1,biny=biny,refy=refy,x.new=x.new,pred.new=pred.new,surv=surv,type=type) #added to the new codel, change the seed to make different results
  temp.1<-NULL
  temp.2<-NULL
  for (l in 1:nx)
  {temp.1<-cbind(temp.1,temp$denm[[l]][,1:ny])
  temp.2<-cbind(temp.2,temp$ie[[l]])}
  te1<-temp$te  #n*nxy
  de1<-temp.1   #n*nxy
  ie1<-temp.2   #n*(nxy*nie)
  cbind(te1,de1,ie1)
  }
  nxy=nx*ny
  te1<-stock(as.matrix(r[,1:nxy]),nrow(x.new))
  de1<-stock(as.matrix(r[,(1+nxy):(2*nxy)]),nrow(x.new))
  nie<-ncol(ie[[1]])
  temp.count<-2*nxy
  for(l in 1:nx)
  {ie2[[l]]<-r[,(temp.count+1):(temp.count+nie)]
  temp.count=temp.count+nie}
  
  if(is.null(w.new))
  {te[2:(n2+1),]<-matrix(apply(te1,2,mean,na.rm=T),ncol=nxy,byrow=T)
  de[2:(n2+1),]<-matrix(apply(de1,2,mean,na.rm=T),ncol=nxy,byrow=T)
  for (l in 1:nx)
    ie[[l]]<-apply(ie2[[l]],2,col_mean1,nrow(x.new))  #first row is the estimated value
  }
  else
  {te[2:(n2+1),]<-matrix(apply(te1,2,weighted.mean,na.rm=T,w=w.new),ncol=nxy,byrow=T)
  de[2:(n2+1),]<-matrix(apply(de1,2,weighted.mean,na.rm=T,w=w.new),ncol=nxy,byrow=T)
  for (l in 1:nx)
    ie[[l]]<-apply(ie2[[l]],2,col_mean1,nrow(x.new),w=w.new)  #first row is the estimated value
  }
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(pred_names,each=ncol(y)),sep=".")
  missing.pred.new<-apply(pred.new,1,anymissing)
  pred.new<-data.frame(pred.new[!missing.pred.new,])
  
  a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=model,
          data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, binpred=F),
          boot.detail=list(pred.new=pred.new,te1=te1,de1=de1,ie1=ie2),w.new=w.new)
  class(a)<-"mma"
  return(a)
}


data<-data.org(x=x,y=y,pred=pred,mediator=mediator,contmed=contmed,binmed=binmed,
               binref=binref,catmed=catmed,catref=catref,jointm=jointm,refy=refy,family1=family1,
               predref=predref,alpha=alpha,alpha2=alpha2,testtype=testtype, w=w)
biny=data$y_type==2
surv=data$y_type==4
if(sum(surv)>0 & is.null(type) & nonlinear)
  type="response"
else if (sum(surv)>0 & is.null(type))
  type="risk"
if(is.null(distn))
  distn<-rep(NA,ncol(data$y))
distn[is.na(distn) & data$y_type==2]="bernoulli"
distn[is.na(distn) & data$y_type==4]="coxph"
distn[is.na(distn) & data$y_type==1]="gaussian"

binpred<-data$binpred
family1<-data$family1

if(binpred) 
{result<-boot.med.binx(data=data,n=n,seed=seed,n2=n2,nonlinear=nonlinear,nu=nu,D=D,distn=distn,family1=family1,
                       w=w,biny=biny,refy=rep(0,ncol(data$y)),surv=surv,type=type)
}
else
{if(is.null(pred.new))
  result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                         D=D,distn=distn,family1=family1,n2=n2,w=w,biny=biny,refy=rep(0,ncol(data$y)),surv=surv,type=type)
else
  result<-boot.med.contx(data=data,margin=margin, n=n,seed=seed, nonlinear=nonlinear,df=df, nu=nu,
                         D=D,distn=distn,family1=family1,n2=n2,w=w,biny=biny,refy=0, x.new=x.new, pred.new=pred.new,surv=surv,type=type,w.new=w.new)
}
result
}



  
