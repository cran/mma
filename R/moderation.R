#Create the moderation function from med, may also be used with mma object, check later
test.moderation<-function(med1,vari,j=1,kx=NULL) #med1 is a med object from the med function, vari is the (vector of) potential moderators 
  #for nonlinear method and the interaction term(s) for linear method
  #j is the jth response if there are multiple responses
  #kx is the kth predictor if k=NULL means all predictor
{test.moderation2<-function(med1,vari,j=1,kx=NULL)
  {binarize<-function(varvec,ref=NULL) #binarize the categorical varvec, ref is the reference group, the first level if null
{a<-factor(varvec)
b=levels(a)
d<-matrix(0,length(varvec),length(b)-1)
if(is.null(ref))
{ref=b[1]
b=b[-1]}
else
  b=b[-grep(ref,b)]
for (k in 1:length(b))
  d[a==b[k],k]<-1
colnames(d)=b
d
}

if(!med1$model$MART)  #if the linear method is used
{temp=Anova(med1$model$model[[j]],type="III")
if(length(vari)==1)
  ln=grep(vari,rownames(temp))
else 
{ln=NULL
for(i in 1:length(vari))
  ln=c(ln,grep(vari[i],rownames(temp)))}
print(temp[ln,])}
else
{temp.name=c(colnames(med1$data$x),colnames(med1$data$dirx))
x=cbind(med1$data$x,med1$data$dirx)
nx=length(temp.name)
#browser()
for(i in 1:length(vari))
{if(is.factor(x[,vari[i]]) | is.character(x[,vari[i]]))
  a=binarize(x[,vari[i]])
else
  a=as.matrix(x[,vari[i]])
if(is.null(kx))
  kx=1:ncol(med1$data$dirx)
for(l in 1:ncol(a))
  for (k in kx)
  {if(is.factor(med1$data$dirx[,k]))
    {temp.dirx=as.numeric(med1$data$dirx[,k])
     temp.dirx=temp.dirx-min(temp.dirx)
     x=cbind(x,a[,l]*temp.dirx)}
   else
    x=cbind(x,a[,l]*med1$data$dirx[,k])
  }
  temp.name=c(temp.name,paste(rep(paste(vari[i],colnames(a),sep=""),each=length(kx)),
                              rep(colnames(med1$data$dirx),ncol(a)),sep="."))
}
colnames(x)=temp.name 
y=med1$data$y[,j]
if(med1$model$Survival[j] & is.null(med1$data$w))
  model=coxph(y~.,data=x)
else if(med1$model$Survival[j])
  model=coxph(y~.,data=x,weights=med1$data$w)
else
  model=glm(y~.,data=x,family=med1$data$family1[[j]],weights=med1$data$w)
temp=Anova(model,type="III")
print(temp[(nx+1):(ncol(x)),])
cat("\nThe H-statistics on MART:\n")
namesdirx=colnames(med1$data$dirx)
for (i in kx)
  for (l in 1:length(vari))
    cat(paste("between ",vari[l]," and ",namesdirx[i],":",sep=""),
        interact.gbm(med1$model$model[[j]],cbind(med1$data$x,med1$data$dirx),i.var=c(namesdirx[i],vari[l])), "\n")
}
}

if(!is.null(med1$a.binx))
{binpred=med1$a.binx$data$binpred
catpred=med1$a.binx$data$catpred
contpred=med1$a.binx$data$contpred
}
else
{binpred=med1$a.contx$data$binpred
  catpred=med1$a.contx$data$catpred
  contpred=med1$a.contx$data$contpred
}

if(is.null(kx))
{if (is.null(binpred))
  for (i in binpred)
    test.moderation2(med1=med1$a.binx,vari=vari,j=j,kx=i)
  if (is.null(contpred))
    for (i in contpred)
      test.moderation2(med1=med1$a.contx,vari=vari,j=j,kx=i)
  if (is.null(catpred))
    for (i in 1:length(catpred))
      test.moderation2(med1=med1$a.binx,vari=vari,j=j,kx=catpred[[i]])
  
}
else{
if(kx%in%binpred)
  test.moderation2(med1=med1$a.binx,vari=vari,j=j,kx=kx)
else if(kx%in%contpred)
  test.moderation2(med1=med1$a.binx,vari=vari,j=j,kx=kx)
else
{z11=rep(F,length(catpred))
for (i in 1:length(catpred))
  z11[i]=kx%in%catpred[[i]]
i=(1:length(catpred))[z11]
 test.moderation2(med1=med1$a.binx,vari=vari,j=j,kx=catpred[[i]])
}}


}
###form the interaction terms
form.interaction<-function(x,pred,inter.cov,predref=NULL,kx=NULL) #create the interaction term.
                                                    #x and binref is the same as in data.org
                                                    #pred is the same set or a subset of pred in data.org, or the mediator vector
                                                    #inter.cov is the name in x that need to form the interaction term
                                                    #kx is the kth predictor if k=NULL means all predictor
{cattobin<-function(x,cat1,cat2=rep(1,length(cat1))) #binaryize the categorical pred
{ ad1<-function(vec)
{vec1<-vec[-1]
vec1[vec[1]]<-1
vec1
}
dim1<-dim(x)
catm<-list(n=length(cat1))
g<-dim1[2]-length(cat1)
ntemp<-colnames(x)[cat1]
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
colnames(f)<-paste(ntemp[j],b,sep="") #changed for error info
hi<-d[d!=l & !is.na(d)]
f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
f[is.na(d),]<-NA
x<-cbind(x,f)
catm<-append(catm,list((g+1):(g+l-1)))
g<-g+length(b)
j<-j+1
}
xname=colnames(x)
x<-x[,-cat1]
x=data.frame(x)
colnames(x)=xname[-cat1]
list(x=x,catm=catm)
}

binarize<-function(varvec,ref=NULL) #binarize the categorical varvec, ref is the reference group, the first level if null
{a<-factor(varvec)
b=levels(a)
if(length(b)==1)
  {if(!is.null(ref))
    if(b==ref)
      return(d=matrix(0,length(varvec),1))
   return(d=matrix(1,length(varvec),1))
  }
d<-matrix(0,length(varvec),length(b)-1)
if(is.null(ref))
{ref=b[1]
b=b[-1]}
else
  b=b[-grep(ref,b)]
for (k in 1:length(b))
  d[a==b[k],k]<-1
colnames(d)=b
d=data.frame(d)
d[is.na(varvec),]=NA
d
}

 pred_names=colnames(pred)
 pred1<-data.frame(pred)
 colnames(pred1)=pred_names
 if(is.null(kx))
   kx=1:ncol(pred1)
 kx1=NULL
 for (i in kx)
  if(nlevels(as.factor(pred1[,i]))==2)
  {if(!is.null(predref))
    pred1[,i]<-ifelse(pred1[,i]==predref,0,1)
   else
   {pred1[,i]<-pred1[,i]
    pred1[,i]<-ifelse(pred1[,i]==levels(pred1[,i])[1],0,1)}
   kx1=c(kx1,i)
  }
 else if(is.character(pred1[,i]) | is.factor(pred1[,i]))
 {pred1[,i]=droplevels(pred1[,i])
  if(!is.null(predref))
   temp.1<-cattobin(data.frame(pred1),i,predref)
  else
   temp.1<-cattobin(data.frame(pred1),i,levels(as.factor(pred1[,i]))[1])
  pred1=temp.1$x
  kx1=c(kx1,temp.1$catm[[2]])
 }
 else
 {kx1=c(kx1,i)}
 temp.name=NULL
 inter=NULL
 #kx=1:ncol(pred1)
 #browser()
 for(i in 1:length(inter.cov))
 {if(is.factor(x[,inter.cov[i]]) | is.character(x[,inter.cov[i]])) #binarize categorical inter.cov
   a=binarize(x[,inter.cov[i]])
  else
   a=as.matrix(x[,inter.cov[i]])
  for(l in 1:ncol(a))
   for (k in kx1)
     #if(is.factor(pred1[,k]))
      # inter=cbind(inter,a[,l]*binarize(pred1[,k])[,1])
     #else
       inter=cbind(inter,a[,l]*pred1[,k])
   temp.name=c(temp.name,paste(rep(paste(inter.cov[i],colnames(a),sep=""),each=length(kx)),
                               rep(colnames(pred1),ncol(a)),sep="."))
   
 }
colnames(inter)=temp.name
inter
}

#estimate and plot the moderate effect from med function
moderate<-function(med1,vari,j=1,kx=1,continuous.resolution=100,plot=T)
{moderate2<-function(med1,vari,j,kx,continuous.resolution,plot)
  {xnames=colnames(med1$data$x)
   pred_names=colnames(med1$data$dirx)
   data1=cbind(med1$data$x,med1$data$dirx)
   colnames(data1)<-c(xnames,pred_names)
   
 if(med1$model$MART)
 {if(is.null(med1$model$type))
   result=plot.gbm(med1$model$model[[j]], i.var=c(pred_names[kx],vari), n.trees=med1$model$best.iter[j],
                   continuous.resolution = continuous.resolution, return.grid=T)
  else
   result=plot.gbm(med1$model$model[[j]], i.var=c(pred_names[kx],vari), n.trees=med1$model$best.iter[j],
                    continuous.resolution = continuous.resolution, return.grid=T,type=med1$model$type)
  if(med1$data$binpred)
  {result=result[(result[,pred_names[kx]]==0 | result[,pred_names[kx]]==1),]
   moderator=unique(result[,2])
   de=matrix(result[,3],2)[2,]-matrix(result[,3],2)[1,]
   result=data.frame(moderator,de)
   if(plot){
     if(is.factor(result$moderator))
       scatterplot(de~moderator, data=result)
     else
       plot(de~moderator, type="l", data=result)}
  }
  else
  {moderator=NULL
   de=NULL
   x=NULL
   for(i in unique(result[,vari]))
   {temp.result=result[result[,vari]==i,]
    x1=temp.result[-nrow(temp.result),1]
    moderator=c(moderator,rep(i,length(x1)))
    de=c(de,diff(temp.result[,3])/diff(temp.result[,1]))
    x=c(x,x1)
   }
   result=data.frame(moderator,de,x)
   if(plot){
     if(is.factor(result$moderator))
       scatterplot(de~x |moderator, smoother=loessLine,data=result)
     else
       levelplot(de~x*moderator, data=result)}
  }
 }
 else
 {model=med1$model$model[[j]]
  pred1=pred_names[kx]
  coef.names=names(model$coefficients)
  beta0=model$coefficients[pred1]  #coefficient for the main dirx
  if(is.na(beta0))
  {pred1=paste(pred1,1,sep="")
   beta0=model$coefficients[pred1]}
  beta=model$coefficients[intersect(grep(pred1,coef.names),grep(vari,coef.names))] #coefficients for the interaction terms
    if(is.factor(med1$data$x[,vari]))
    result=data.frame(moderator=c("ref",names(beta)),de=c(beta0,beta0+beta))
  else if (med1$data$binpred){
    if(length(beta)==1)
     result1=data.frame(moderator=sort(unique(med1$data$x[,vari])),de=beta0+beta*sort(unique(med1$data$x[,vari])))
    else
    {temp.order=order(c(0,med1$data$x[med1$data$dirx[,kx]==1,vari]))
     de=c(beta0,beta0+as.matrix(med1$data$x[med1$data$dirx[,kx]==1,intersect(grep(pred1,xnames),grep(vari,xnames))])%*%beta)
     result1=data.frame(moderator=c(0,med1$data$x[med1$data$dirx[,kx]==1,vari])[temp.order],de=de[temp.order])
    }
     temp.2=NULL
     for (i in unique(result1$moderator))
       temp.2=c(temp.2,result1[result1$moderator==i,"de"][1])
    result=data.frame(moderator=unique(result1$moderator),de=temp.2)
  }
  else{
    if(length(beta)==1)
      result=data.frame(moderator=sort(unique(med1$data$x[,vari])),de=beta0+beta*sort(unique(med1$data$x[,vari])))
    else
    {temp.in=!is.na(med1$data$dirx[,kx]) & med1$data$dirx[,kx]!=0 & !is.na(med1$data$x[,vari])
     temp.order=order(med1$data$x[temp.in,vari])
     de=beta0+diag(1/med1$data$dirx[temp.in,kx])%*%as.matrix(med1$data$x[temp.in,intersect(grep(pred1,xnames),grep(vari,xnames))])%*%beta
     result1=data.frame(moderator=med1$data$x[temp.in,vari][temp.order],de=de[temp.order])
     temp.2=NULL
     for (i in unique(result1$moderator))
       temp.2=c(temp.2,result1[result1$moderator==i,"de"][1])
     result=data.frame(moderator=unique(result1$moderator),de=temp.2)}
  }
  if(plot){
    if(is.factor(result$moderator))
      scatterplot(de~moderator, data=result)
    else
      plot(de~moderator, type="l", data=result)}
}
a=list(result=result,med1=med1,vari=vari,j=j,kx=kx)
class(a)="moderate"
a
}
if(!is.null(med1$a.binx))
{binpred=med1$a.binx$binpred
 catpred=med1$a.binx$catpred
 contpred=med1$a.binx$contpred
}
else
{binpred=med1$a.contx$binpred
 catpred=med1$a.contx$catpred
 contpred=med1$a.contx$contpred
}

if(kx %in% contpred)
  moderate2(med1=med1$a.contx,vari=vari,j=j,kx=kx,
            continuous.resolution=continuous.resolution,plot=plot)
else
  moderate2(med1=med1$a.binx,vari=vari,j=j,kx=kx,
            continuous.resolution=continuous.resolution,plot=plot)
}

#make inferences on moderation (mediated or not) effects from the mma function.
boot.mod<-function(mma1,vari,continuous.resolution=10,
                   w=NULL,n=20,
                   x.new=NULL,w.new=NULL,pred.new=NULL,cova.new=NULL,xj=1,margin=1,xmod=vari,df1=1)
  #boots=T for bootstrap method
  #continuous.resolution: for continuous moderator, this is the number of points to be taken from 
  ##min to max by 1/continuous.resolution. For categorical moderator, this is the categories to moderate, 
  ##all if it is not set. If there is no enough case with the 1/continuous.resolution quintile, error shows
  ##to reduce continuous.resolution.
  #kx and jy can be vectors #kx should be xj
{boot.mod.binx<-function(mma1,vari,plot=T,continuous.resolution=100,n2=NULL,
                         n=20,w=rep(1,nrow(mma1$data$x)),xj=1,xmod=vari)
  #n2 is the time of bootstrap if set as null. It has to be less or equal to the number of bootstrap
{mod.binx<-function(vari,continuous.resolution,n,x,y,dirx,contm,catm,
                    jointm,cova,allm,full.model,best.iter1,surv,type,w=w,moder.level1=NULL,xj=1,xmod) #
{xnames<-colnames(x)
pred_names<-colnames(dirx)  

te.binx<-function(full.model,new1,new0,best.iter1=NULL,surv,type)       
{te<-NULL
for(m in 1:length(full.model))
  if(surv[m] & !is.null(best.iter1[m]))
    te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
  else if (surv[m])
    te[m]<-mean(predict(full.model[[m]],new1,type=type),na.rm=T)- mean(predict(full.model[[m]],new0,type=type),na.rm=T)
  else if(!is.null(best.iter1[m]))
    te[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
  else
    te[m]<-mean(predict(full.model[[m]],new1),na.rm=T)- mean(predict(full.model[[m]],new0),na.rm=T)
  te
}

med.binx.contm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,xmod,xnames)  
{n3<-nrow(nom1)+nrow(nom0)
marg.m<-c(nom1[,med],nom0[,med])[sample(1:n3,replace=T)]
new1<-nom1
new1[,med]<-marg.m[1:nrow(nom1)]
new0<-nom0
new0[,med]<-marg.m[(nrow(nom1)+1):n3]
if(!is.null(xmod))
{temp.x=intersect(grep(xnames[med],xnames),grep(xmod,xnames))
if(sum(temp.x)>0)
{m.t=1
m.t2=form.interaction(new0,new0[,med],inter.cov=xmod)
m.t3=form.interaction(new1,new1[,med],inter.cov=xmod)
for (m.t1 in temp.x)
{new0[,m.t1]=m.t2[,m.t]
new1[,m.t1]=m.t3[,m.t]
m.t=m.t+1}}
}
dir.nom<-NULL
for(m in 1:length(full.model))
  if(surv[m] & !is.null(best.iter1[m]))
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
else if(surv[m])
  dir.nom[m]<-mean(predict(full.model[[m]],new1,type=type),na.rm=T)- mean(predict(full.model[[m]],new0,type=type),na.rm=T)
else if(!is.null(best.iter1[m]))
  dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
else
  dir.nom[m]<-mean(predict(full.model[[m]],new1),na.rm=T)- mean(predict(full.model[[m]],new0),na.rm=T)
dir.nom
}

med.binx.jointm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,temp.rand,xmod,xnames)  
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
  if(!is.null(xmod))
    for (z in med)
    {temp.x=intersect(grep(xnames[z],xnames),grep(xmod,xnames))
    if(sum(temp.x)>0)
    {m.t=1
    m.t2=form.interaction(new0,new0[,z],inter.cov=xmod)
    m.t3=form.interaction(new1,new1[,z],inter.cov=xmod)
    for (m.t1 in temp.x)
    {new0[,m.t1]=m.t2[,m.t]
    new1[,m.t1]=m.t3[,m.t]
    m.t=m.t+1}}
    }
  dir.nom<-NULL
  for (m in 1:length(full.model))
    if(surv[m] & !is.null(best.iter1[m]))
      dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T)
  else if(surv[m])
    dir.nom[m]<-mean(predict(full.model[[m]],new1,type=type),na.rm=T)- mean(predict(full.model[[m]],new0,type=type),na.rm=T)
  else if(!is.null(best.iter1[m]))
    dir.nom[m]<-mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T)
  else
    dir.nom[m]<-mean(predict(full.model[[m]],new1),na.rm=T)- mean(predict(full.model[[m]],new0),na.rm=T)
  dir.nom
}

med.binx.catm<-function(full.model,nom1,nom0,med,best.iter1=NULL,surv,type,xmod,xnames)  
{n3<-nrow(nom1)+nrow(nom0)
temp.rand<-unlist(list(nom1[,med],nom0[,med]))[sample(1:n3,replace=T)]
marg.m1<-temp.rand[1:nrow(nom1)]
marg.m2<-temp.rand[(nrow(nom1)+1):n3]
dir.nom<-rep(0,length(full.model))
for (m in 1:length(full.model))
  for (i in levels(x[,med]))
  {new1<-nom1
  new1[1:dim(new1)[1],med]<-i
  new0<-nom0
  new0[1:dim(new0)[1],med]<-i
  if(!is.null(xmod))
  {temp.x=intersect(grep(xnames[med],xnames),grep(xmod,xnames))
  if(sum(temp.x)>0)
  {m.t=1
  m.t2=form.interaction(new0,new0[,med],inter.cov=xmod)
  m.t3=form.interaction(new1,new1[,med],inter.cov=xmod)
  for (m.t1 in temp.x)
  {new0[,m.t1]=m.t2[,m.t]
  new1[,m.t1]=m.t3[,m.t]
  m.t=m.t+1}}
  }
  p<-mean(temp.rand==i,na.rm=T)
  if(surv[m] & !is.null(best.iter1[m]))
    dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m],type=type),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m],type=type),na.rm=T))
  else if(surv[m])
    dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,type=type),na.rm=T)- mean(predict(full.model[[m]],new0,type=type),na.rm=T))
  else if(!is.null(best.iter1[m]))
    dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1,best.iter1[m]),na.rm=T)- mean(predict(full.model[[m]],new0,best.iter1[m]),na.rm=T))
  else
    dir.nom[m]<-dir.nom[m]+p*(mean(predict(full.model[[m]],new1),na.rm=T)- mean(predict(full.model[[m]],new0),na.rm=T))
  }
dir.nom
}

#1.get the model
x2<-cbind(x,dirx)
colnames(x2)<-c(xnames,pred_names)

#1.5 prepare for the moderator
if(is.null(moder.level1)){
  moder.level=NULL
  if(is.factor(x[,vari]))
  {if(continuous.resolution==100)
    moder.level=levels(x[,vari])
  else
    moder.level=continuous.resolution
  for (i in moder.level)
  {temp.all=(data$x[,vari]==i)
  if(sum(apply(as.matrix(dirx[temp.all,]==1),2,sum,na.rm=T)==0)>1 | 
     sum(dirx[temp.all,]==1,na.rm=T)==length(dirx[temp.all,1][!is.na(dirx[temp.all,1])]))
    stop("Error: need to reduce the continuous.resolution") #error if the group has all dirx=0 or 1
  }
  temp.q=NULL
  }
  else
  {temp.q=quantile(unique(x[,vari]),probs=(seq(0,1,by=1/continuous.resolution))[-1],na.rm=T)  #add unique to take care of repeats
  for(i in 1:length(temp.q))
  {if (i==1)
    temp.all=(x[,vari]<=temp.q[i])
  else
    temp.all=(x[,vari]<=temp.q[i] & x[,vari]>temp.q[i-1])
if(sum(apply(as.matrix(dirx[temp.all,]==0),2,sum,na.rm=T)==0)>1 | 
     sum(dirx[temp.all,]==0,na.rm=T)==length(dirx[temp.all,1][!is.na(dirx[temp.all,1])]))
    stop("Error: need to reduce the continuous.resolution") #error if the group has all dirx=0 or 1
  #if(!is.null(w))
  #{w.moder=c(w.moder,sum(w[temp.all]))
  # moder.level=c(moder.level, weighted.mean(data$x[temp.all,vari],w[temp.all]))}
  #else}
  moder.level=c(moder.level,mean(data$x[temp.all,vari],na.rm=T))
  }}}
else
  {moder.level=moder.level1$moder.level
   temp.q=moder.level1$cont.moder.q}

nmod=length(moder.level)
#2. prepare for the store of results
#set.seed(seed)
te<-matrix(NA,n,ncol(y)*nmod)
colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(moder.level,each=ncol(y)),sep=".")
if(!is.null(jointm))
{denm<-matrix(NA,n,ncol(y)*(1+length(c(contm,catm))+jointm[[1]]))
dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",colnames(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")
}
else
{denm<-matrix(NA,n,ncol(y)*(1+length(c(contm,catm))))
dimnames(denm)[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",colnames(x)[c(contm,catm)]),each=ncol(y)),sep=".")
}
denm<-rep(list(denm),nmod)
ie<-denm

#3. repeat to get the mediation effect
for(q1 in 1:length(moder.level)){
  if(is.factor(x[,vari]))
    temp.all=(x[,vari]==moder.level[q1] & !is.na(x[,vari]))
  else if (q1==1)
    temp.all=(x[,vari]<=temp.q[q1] & !is.na(x[,vari]))
  else
    temp.all=(x[,vari]<=temp.q[q1] & x[,vari]>temp.q[q1-1] & !is.na(x[,vari]))
  
  dirx1=data.frame(dirx[temp.all,])
  names(dirx1)=pred_names
  w.temp=w[temp.all]
  x2.1<-x2[temp.all,]
  x2.2=data.frame(x[temp.all,])
  colnames(x2.2)=xnames
  colnames(x2.1)<-c(xnames,pred_names)
 # 
  for (k in 1:n)
  {#3.1 get the te         full.model,x,y,dirx,best.iter1=NULL
    x0.temp<-apply(as.matrix(dirx1[,xj]==1),1,sum)==0  #indicator of the reference group
    if(sum(x0.temp)==0) break #to break out if there is not reference group
    x0<-x2.1[x0.temp,]
    if(is.null(w.temp))
    {w1<-NULL
    w0<-NULL}
    else
      w0<-w.temp[x0.temp]
    
    for (l in xj)  #l indicate the lth predictor
    {if(sum(dirx1[,l]==1)==0) next  #next if there is no this group
      x1.2<-x2.1[dirx1[,l]==1,]
      if(!is.null(w.temp))
        w1<-w.temp[dirx1[,l]==1]
      new1<-x1.2[sample(1:nrow(x1.2),replace=T,prob=w1),] #floor(n3/2),
      new0<-x0[sample(1:nrow(x0),replace=T,prob=w0),] #floor(n3/2),
      
      if(!is.null(xmod)  & !is.factor(x[,xmod]))
        for(z in allm){
          temp.x=intersect(grep(xnames[z],xnames),grep(xmod,xnames))
          if(sum(temp.x)>0)
          {m.t=1
          m.t2=form.interaction(new0,new0[,z],inter.cov=xmod)
          m.t3=form.interaction(new1,new1[,z],inter.cov=xmod)
          for (m.t1 in temp.x)
          {new0[,m.t1]=m.t2[,m.t]
          new1[,m.t1]=m.t3[,m.t]
          m.t=m.t+1}}
        }
      
      te[k,((q1-1)*ncol(y)+1):(q1*ncol(y))]<-te.binx(full.model,new1,new0,best.iter1,surv,type)  
      temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
      #the indirect effect of all mediators
      temp.ie<-te[k,((q1-1)*ncol(y)+1):(q1*ncol(y))]-med.binx.jointm(full.model,new1,new0,allm,best.iter1,surv,type,temp.rand,xmod,xnames) #add temp.rand
      #new method to calculate the direct effect     
      x.temp=rbind(x2.2[dirx1[,l]==1,],x2.2[x0.temp,])
      new1.temp=cbind(x.temp[temp.rand[1:nrow(x1.2)],],dirx1[dirx1[,l]==1,])
      new0.temp=cbind(x.temp[temp.rand[(nrow(x1.2)+1):(nrow(x1.2)+nrow(x0))],],dirx1[x0.temp,])
      colnames(new1.temp)<-c(xnames,pred_names)
      colnames(new0.temp)<-c(xnames,pred_names)
      if(!is.null(xmod) & !is.factor(x[,xmod])){
        temp.x=intersect(grep(pred_names[l],xnames),grep(xmod,xnames))
        if(sum(temp.x)>0)
        {m.t=1
        m.t2=form.interaction(new0.temp,dirx1[x0.temp,],inter.cov=xmod)
        m.t3=form.interaction(new1.temp,dirx1[dirx1[,l]==1,],inter.cov=xmod)

        for (m.t1 in temp.x)
        {new0.temp[,m.t1]=m.t2[,m.t]
        new1.temp[,m.t1]=m.t3[,m.t]
        m.t=m.t+1}}}
      denm[[q1]][k,1:ncol(y)]<-te.binx(full.model,new1.temp,new0.temp,best.iter1,surv,type) #add temp.rand
      
      j<-2
      #3.2 mediation effect from the continuous mediator
      if (!is.null(contm))
        for (i in contm)          #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[q1]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.contm(full.model,new1,new0,i,best.iter1,surv,type,xmod,xnames)
        j<-j+1}
      #3.3.mediation effect from the categorical mediator
      if (!is.null(catm))
        for (i in catm)           #full.model,x,y,med,dirx,best.iter1=NULL
        {denm[[q1]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.catm(full.model,new1,new0,i,best.iter1,surv,type,xmod,xnames)
        j<-j+1}
      #3.4 mediation effect from the joint mediators
      if (!is.null(jointm))
        for (i in 1:jointm[[1]])          #full.model,x,y,med,dirx,best.iter1=NULL
        {temp.rand<-sample(1:(nrow(x1.2)+nrow(x0)),replace=T)# no need for:prob=c(w1,w0) --redundant
        denm[[q1]][k,(ncol(y)*(j-1)+1):(ncol(y)*j)]<-med.binx.jointm(full.model,new1,new0,jointm[[i+1]],best.iter1,surv,type,temp.rand,xmod,xnames)
        j<-j+1}
      #3.5 recalculate the total effect and get the indirect effects
      ie[[q1]][k,]<-te[k,((q1-1)*ncol(y)+1):(q1*ncol(y))]-denm[[q1]][k,]
      ie[[q1]][k,1:ncol(y)]<-temp.ie
      te[k,((q1-1)*ncol(y)+1):(q1*ncol(y))]<-denm[[q1]][k,1:ncol(y)]+temp.ie
      if(!is.null(jointm))
        dimnames(ie[[q1]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",colnames(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ncol(y)),sep=".")#c("all",colnames(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep=""))
      else
        dimnames(ie[[q1]])[[2]]<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",colnames(x)[c(contm,catm)]),each=ncol(y)),sep=".") #c("all",colnames(x)[c(contm,catm)])
    }
  }}
names(denm)<-moder.level
names(ie)<-moder.level
a<-list(denm=denm,ie=ie,te=te,moder.level=list(moder.level=moder.level,cont.moder.q=temp.q,moder=x[,vari]),data=data,mod=T)
class(a)<-"med"
return(a)
}
#browser()
data=mma1$data
x=data$x
y=data$y
dirx=data$dirx
contm=data$contm
catm=data$catm
jointm=data$jointm
cova=data$cova
allm=c(contm,catm)
if (is.null(allm))
  stop("Error: no potential mediator is specified")
xnames<-colnames(x)
pred_names<-colnames(dirx)
ynames=colnames(y)

full.model<-mma1$model$model
best.iter1<-mma1$model$best.iter
surv<-mma1$model$Survival
type<-mma1$model$type

temp<-mod.binx(vari=vari,continuous.resolution=continuous.resolution,n=n,x=x,y=y,dirx=dirx,
               contm=contm,catm=catm,jointm=jointm,cova=cova,allm=allm,full.model=full.model,
               best.iter1=best.iter1,surv=surv,type=type,w=w,moder.level1=NULL,xj=xj,xmod=xmod)

ny=ncol(y)
nx=1
nmod=length(temp$moder.level$moder.level)
if(is.null(n2))
  n2=ifelse(is.null(mma1$all_boot),0,nrow(mma1$all_boot))
te<-matrix(0,n2+1,ny*nmod)
de<-matrix(0,n2+1,ny*nmod)
if(is.null(jointm))
{ie<-matrix(0,n2,ny*(1+length(c(contm,catm))))
ie1<-matrix(0,nmod,ny*(1+length(c(contm,catm))))
colnames(ie)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[c(contm,catm)]),each=ny),sep=".")
#dimnames(ie)[[3]]=dimnames(temp$te)[[3]]
colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[c(contm,catm)]),each=ny),sep=".")
rownames(ie1)<-temp$moder.level$moder.level
#dimnames(ie1)[[3]]=dimnames(temp$te)[[3]]
}
else 
{ie<-matrix(0,n2,ny*(1+length(c(contm,catm))+jointm[[1]]))
colnames(ie)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
ie1<-matrix(0,nmod,ny*(1+length(c(contm,catm))+jointm[[1]]))
colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[c(contm,catm)],paste("j",1:jointm[[1]],sep="")),each=ny),sep=".")
rownames(ie1)<-temp$moder.level$moder.level
#dimnames(ie1)[[3]]=dimnames(temp$te)[[3]]
}
ie<-rep(list(ie),nmod)
names(ie)<-temp$moder.level$moder.level

te[1,]<-apply(temp$te,2,mean,na.rm=T)
temp.1<-temp$te
for (l in 1:nmod)
{temp.1[,l]<-temp$denm[[l]][,1:ny]
ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)}  #first row is the estimated value
de[1,]<-apply(temp.1,2,mean,na.rm=T)

moder.level1=temp$moder.level

if(n2>0){
  for (i in 1:n2)
  {boots<-mma1$all_boot[i,]
  x1<-data.frame(x[boots,])
  y1<-data.frame(y[boots,])
  wz=w[boots]
  pred1<-data.frame(dirx[boots,])
  full.model=mma1$all_model[[i]]
  best.iter1=mma1$all_iter[i,]
  colnames(x1)=xnames
  colnames(y1)=ynames
  colnames(pred1)=pred_names
  
  temp<-mod.binx(vari,continuous.resolution,n,x1,y1,pred1,contm,catm,
                 jointm,cova,allm,full.model,best.iter1,surv,type,wz,moder.level1,xj,xmod)
  
  te[1+i,]<-apply(temp$te,2,mean,na.rm=T)
  temp.1<-temp$te
  for (l in 1:nmod)
  {temp.1[,l]<-temp$denm[[l]][,1:ny]
  ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)}  #first row is the estimated value
  de[1+i,]<-apply(temp.1,2,mean,na.rm=T)
  print(i)
  }}

moder.level=moder.level1$moder.level

te1=matrix(te[-1,],n2,ny*nmod)
de1=matrix(de[-1,],n2,ny*nmod)
colnames(te1)<-paste(paste("y",1:ncol(y),sep=""),rep(moder.level,each=ncol(y)),sep=".")
colnames(de1)<-paste(paste("y",1:ncol(y),sep=""),rep(moder.level,each=ncol(y)),sep=".")
colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(moder.level,each=ncol(y)),sep=".")
colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(moder.level,each=ncol(y)),sep=".")

a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te1,de=de1), 
        data=list(x=x,y=y,dirx=dirx,contm=contm,catm=catm,jointm=jointm,binpred=T),model=mma1$model,
        moder.level=moder.level1,mod=T)
class(a)<-"mma"
return(a)
}

boot.mod.contx<-function(mma1,vari,continuous.resolution=10,
                         w=rep(1,nrow(mma1$data$x)),n=20,
                         x.new=NULL,w.new=NULL,
                         pred.new=NULL,cova.new=NULL,xj=1,df1=1,xmod=vari,margin=1)
{mod.contx<-function(vari,continuous.resolution,x,y,dirx,binm,contm,catm,jointm,cova, n,x.new=x,
                     pred.new=dirx, cova.new=cova, w=rep(1,nrow(x)), w.new=w,full.model,best.iter1,
                     surv,type,moder.level1,nonlinear=nonlinear,df1=1,n2=NULL,xj=1,xmod=vari,margin=1)
  
{if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")
  
  xnames<-colnames(x)
  pred_names<-colnames(dirx)
  ynames<-colnames(y)
  cova_names<-colnames(cova)
  
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
  
  
  dist.m.given.x<-function(x,dirx,binm=NULL,contm=NULL,catm=NULL,nonlinear,df1,w,cova) #give the model and residual of m given x
  {
    getform=function(z,nonlinear,df1)
    {if(!nonlinear)
      formu="x[,i]~."
    else
    {names.z=colnames(z)
    temp.t=unlist(lapply(z,is.character)) | unlist(lapply(z,is.factor))
    names.z1=names.z[!temp.t]
    names.z2=names.z[temp.t]
    if(length(names.z1)==0)
      formu="x[,i]~."
    else if (length(names.z2)==0)
      formu=paste("x[,i]~",paste(paste("ns(",names.z1,",","df=",df1,")",sep=""),collapse="+"),sep="")
    else
      formu=paste("x[,i]~",paste(paste("ns(",names.z1,",","df=",df1,")",sep=""),collapse="+"),"+",
                  paste(names.z2,collapse="+"),sep="")
    }
    formu
    }
    #browser()  
    models<-NULL
    x=data.frame(x)
    res<-NULL
    temp.namec=colnames(x)
    indi=NULL                               #indi indicate if not all mediators, the columns of mediators that needs covariates
    if(!is.null(cova))
      if(length(grep("for.m",names(cova)))!=0)
        for (i in 1:length(cova[[2]]))
          indi=c(indi,grep(cova[[2]][i],temp.namec))
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm<-c(binm,catm[[i]])}
    
    z<-dirx
    z.name=paste("predictor",1:ncol(z),sep=".")
    colnames(z)=z.name
    
    if(!is.null(cova))
    {if (length(grep("for.m",names(cova)))==0)#create the predictor matrix z
      z<-cbind(z,cova)
    else if(length(grep("for.m",names(cova)))!=0)
    {
      z1<-cbind(z,cova[[1]])
      form1=getform(z1,nonlinear,df1)
    }}
    
    form0=getform(z,nonlinear,df1)
    j<-1
    
    if(!is.null(binm))
    {for(i in binm)
    {if(!i%in%indi)
    {models[[j]]<-glm(form0,data=data.frame(z),family=binomial(link = "logit"),weights=w)
    res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))}
      else
      {models[[j]]<-glm(form1,data=data.frame(z1),family=binomial(link = "logit"),weights=w)
      res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z1)))}
      j<-j+1}
    }
    for (i in contm)
    {if(!i%in%indi)
      models[[j]]<-glm(as.formula(form0),data=data.frame(z),family=gaussian(link="identity"),weights=w)
    else
      models[[j]]<-glm(as.formula(form1),data=data.frame(z1),family=gaussian(link="identity"),weights=w)
    res<-cbind(res,models[[j]]$res)
    j<-j+1
    }
    list(models=models,varmat=var(res))
  }
  
  
  sim.xm<-function(distmgivenx,x1,dirx,binm,contm,catm,nonlinear,df1,cova)  #added nonlinear and df1 to sim.xm
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
  
  x1=data.frame(x1)
  temp.namec=colnames(x1)
  indi=NULL                               #indi indicate if not all mediators, the columns of mediators that needs covariates
  if(!is.null(cova))
    if(length(grep("for.m",names(cova)))!=0)
     for (i in 1:length(cova[[2]]))
      indi=c(indi,grep(cova[[2]][i],temp.namec))
  
  means<-NULL
  z<-dirx
  z.name=paste("predictor",1:ncol(z),sep=".")
  colnames(z)=z.name
  
  if(!is.null(cova))
  {if(length(grep("for.m",names(cova)))==0)   #create the predictor matrix z
    z<-cbind(z,cova)
  else 
    z1<-cbind(z,cova[[1]])}
  
  binm1<-binm
  
  if(!is.null(catm))
  {for (i in 2:(catm$n+1))
    binm1<-c(binm1,catm[[i]])}
  if(!is.null(binm1))
    for (i in 1:length(binm1))
    {if(binm1[i]%in%indi)
      means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z1)))
    else  
      means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z)))}
  if(!is.null(contm))
    for (i in (length(binm1)+1):length(c(binm1,contm)))
    {if(contm[i-length(binm1)]%in%indi)
      means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z1)))
    else
      means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z)))}
  
  if(dim(means)[2]==1)                                                   #added in the new program, in case there is only one mediator
  {sim.m<-suppressWarnings(rnorm(length(means),mean=means,sd=sqrt(distmgivenx$varmat)))     #added in the new program
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
  x<-data.frame(x[nonmissing,])
  colnames(x)=xnames
  y<-data.frame(y[nonmissing,])
  if(!is.null(cova))
  {cova=data.frame(cova[nonmissing,])
  colnames(cova)=cova_names}
  colnames(y)<-ynames
  pred<-data.frame(dirx[nonmissing,])
  colnames(pred)<-pred_names
  w<-w[nonmissing]

  
  #2. prepare for the store of results
  #set.seed(seed)
  #n.new1=sum(mod.level1$n.moder.level)-length(mod.level1$n.moder.level)
  #  te<-matrix(0,n.new,ncol(dirx)*ncol(y))
  
  #3. get the joint distribution of m given x
 # browser()
  distmgivenx<-dist.m.given.x(x,pred,binm,contm,catm,nonlinear,df1,w,cova)
  te1.0<-NULL
  denm1.0<-NULL
  denm1.1<-NULL
  
  n1<-dim(x)[1]
  nmod=moder.level1$n.moder.level
  
  #4. repeat to get the mediation effect
  for (l in 1:nmod) {    #browser()
    level=moder.level1$levels[,l]
    x.new1=data.frame(x.new[level,])
    colnames(x.new1)=xnames
    n.new=nrow(x.new1)
    pred.new1=data.frame(pred.new[level,])
    colnames(pred.new1)=pred_names
    if(!is.null(cova.new))
    {cova.new1=data.frame(cova.new[level,])
     colnames(cova.new1)=cova_names}
    else
      cova.new1=NULL
    denm1<-NULL
    denm1.2=NULL
    te1<-NULL
    for (k in 1:n)
    {new0<-sim.xm(distmgivenx,x.new1,pred.new1,binm,contm,catm,nonlinear,df1,cova.new1) #draw ms conditional on x.new
    temp.pred<-pred.new1
    temp.pred[,xj]<-temp.pred[,xj]+margin
    if(!is.null(xmod))   #allows the interaction of pred with xmod
    {cova.new2=cova.new1
    x.new2=x.new1
    if(!is.null(cova.new1))
    {temp.cova=intersect(grep(pred_names[xj],cova_names),grep(xmod,cova_names))
    if(sum(temp.cova)>0)
    {m.t=1
    #browser()
    m.t2=form.interaction(cova.new1,temp.pred[,xj],inter.cov=xmod)
    for (m.t1 in temp.cova)
    {cova.new2[,m.t1]=m.t2[,m.t]
    m.t=m.t+1}
    }}
    temp.x=intersect(grep(pred_names[xj],xnames),grep(xmod,xnames))
    if(sum(temp.x)>0)
    {m.t=1
    m.t2=form.interaction(x.new1,temp.pred[,xj],inter.cov=xmod)
    for (m.t1 in temp.x)
    {x.new2[,m.t1]=m.t2[,m.t]
    m.t=m.t+1}}
    new1<-sim.xm(distmgivenx,x.new2,temp.pred,binm,contm,catm,nonlinear,df1,cova.new2)  #draw from the conditional distribution of m given x
    }
    else
      new1<-sim.xm(distmgivenx,x.new1,temp.pred,binm,contm,catm,nonlinear,df1,cova.new1)  #draw from the conditional distribution of m given x
    new1<-cbind(new1,temp.pred)   #draw ms conditional on x.new+margin
    new0<-cbind(new0,pred.new1)
    
    if(!is.null(xmod))
      for(z in allm){
        temp.x=intersect(grep(xnames[z],xnames),grep(xmod,xnames))
        if(sum(temp.x)>0)
        {m.t=1
        m.t2=form.interaction(new0,new0[,z],inter.cov=xmod)
        m.t3=form.interaction(new1,new1[,z],inter.cov=xmod)
        for (m.t1 in temp.x)
        {new0[,m.t1]=m.t2[,m.t]
        new1[,m.t1]=m.t3[,m.t]
        m.t=m.t+1}}
      }

    denm2<-NULL
    
   # browser()
    
    sample.temp<-sample(1:n.new,2*n.new,replace = T,prob=w.new[level])   #random sample from the original data
    
    #4.0.0 get the total indirect effect
    temp.new1<-new1
    temp.new1[,allm]<-x.new1[sample.temp[1:n.new],allm]
    temp.new0<-new0
    temp.new0[,allm]<-x.new1[sample.temp[(n.new+1):(2*n.new)],allm]
    if(!is.null(xmod))
      for(z in allm){
        temp.x=intersect(grep(xnames[z],xnames),grep(xmod,xnames))
        if(sum(temp.x)>0)
        {m.t=1
        m.t2=form.interaction(x.new1[sample.temp[1:n.new],],x.new1[sample.temp[1:n.new],z],inter.cov=xmod)
        m.t3=form.interaction(x.new1[sample.temp[(n.new+1):(2*n.new)],],x.new1[sample.temp[(n.new+1):(2*n.new)],z],inter.cov=xmod)
        for (m.t1 in temp.x)
        {temp.new1[,m.t1]=m.t2[,m.t]
        temp.new0[,m.t1]=m.t3[,m.t]
        m.t=m.t+1}}
      }
    
    for (m in 1:ncol(y))
      if(surv[m] & !is.null(best.iter1[m]))
        denm3<-(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin
    else if(surv[m])
      denm3<-(predict(full.model[[m]],temp.new1,type=type)-predict(full.model[[m]],temp.new0,type=type))/margin
    else
      denm3<-(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin
    
    #4.0 get the direct effect
    temp.new1<-x.new[sample.temp[1:n.new],]
    temp.new1=cbind(temp.new1,temp.pred)
    temp.new0<-x.new[sample.temp[(n.new+1):(2*n.new)],]
    temp.new0=cbind(temp.new0,pred.new1)
    colnames(temp.new1)<-c(xnames,pred_names)
    colnames(temp.new0)<-c(xnames,pred_names)
    
    if(!is.null(xmod)){
      temp.x=intersect(grep(pred_names[l],xnames),grep(xmod,xnames))
      if(sum(temp.x)>0)
      {m.t=1
      m.t2=form.interaction(temp.new1,temp.pred[,l],inter.cov=xmod)
      m.t3=form.interaction(temp.new0,pred.new1[,l],inter.cov=xmod)
      for (m.t1 in temp.x)
      {temp.new1[,m.t1]=m.t2[,m.t]
      temp.new0[,m.t1]=m.t3[,m.t]
      m.t=m.t+1}}
    }
    
    for (m in 1:ncol(y))
      if(surv[m] & !is.null(best.iter1[m]))
        denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m],type=type)-predict(full.model[[m]],temp.new0,best.iter1[m],type=type))/margin)
    else if(surv[m])
      denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,type=type)-predict(full.model[[m]],temp.new0,type=type))/margin)
    else
      denm2<-cbind(denm2,(predict(full.model[[m]],temp.new1,best.iter1[m])-predict(full.model[[m]],temp.new0,best.iter1[m]))/margin)
    
    #4.1 get the te
    te0<-NULL
    for(m in 1:ncol(y))
      if(surv[m] & !is.null(best.iter1[m]))
        te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m],type=type)-predict(full.model[[m]],new0,best.iter1[m],type=type))/margin)
    else if(surv[m])
      te0<-c(te0, (predict(full.model[[m]],new1,type=type)-predict(full.model[[m]],new0,type=type))/margin)
    else
      te0<-c(te0, (predict(full.model[[m]],new1,best.iter1[m])-predict(full.model[[m]],new0,best.iter1[m]))/margin)
    te1<-cbind(te1,te0)
    
    #4.2 mediation effect from the single mediator
    # browser()
    if (!is.null(listm$single))
      for (i in 1:length(listm$single))
      {new1.nm<-new1
      new0.nm<-new0
      temp.m<-x.new1[sample.temp,listm$single[i]]
      new1.nm[,listm$single[i]]<-temp.m[1:n.new]    #draw m from its original distribution
      new0.nm[,listm$single[i]]<-temp.m[(n.new+1):(2*n.new)]    #draw m from its original distribution
      
      if(!is.null(xmod))
      {temp.x=intersect(grep(xnames[listm$single[i]],xnames),grep(xmod,xnames))
      if(sum(temp.x)>0)
      {m.t=1
      m.t2=form.interaction(new1.nm,new1.nm[,listm$single[i]],inter.cov=xmod)
      m.t3=form.interaction(new0.nm,new0.nm[,listm$single[i]],inter.cov=xmod)
      for (m.t1 in temp.x)
      {new1.nm[,m.t1]=m.t2[,m.t]
      new0.nm[,m.t1]=m.t3[,m.t]
      m.t=m.t+1}}
      }
      
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,type=type)-predict(full.model[[m]],new0.nm,type=type))/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
      }
    
    #4.3.mediation effect from the joint mediator
    if (!is.null(listm$multi))
      for (i in 2:(listm$multi[[1]]+1))
      {new1.nm<-new1
      new0.nm<-new0
      new1.nm[,listm$multi[[i]]]<-x.new1[sample.temp[1:n.new],listm$multi[[i]]]    #draw m from its original distribution
      new0.nm[,listm$multi[[i]]]<-x.new1[sample.temp[(n.new+1):(2*n.new)],listm$multi[[i]]]    #draw m from its original distribution
      
      if(!is.null(xmod))
        for (z in listm$multi[[i]])
        {temp.x=intersect(grep(xnames[z],xnames),grep(xmod,xnames))
        if(sum(temp.x)>0)
        {m.t=1
        m.t2=form.interaction(new1.nm,new1.nm[,z],inter.cov=xmod)
        m.t3=form.interaction(new0.nm,new0.nm[,z],inter.cov=xmod)
        for (m.t1 in temp.x)
        {new1.nm[,m.t1]=m.t2[,m.t]
        new0.nm[,m.t1]=m.t3[,m.t]
        m.t=m.t+1}}
        }
      
      for(m in 1:ncol(y))
        if(surv[m] & !is.null(best.iter1[m]))
          denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m],type=type)-predict(full.model[[m]],new0.nm,best.iter1[m],type=type))/margin)
      else if(surv[m])
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,type=type)-predict(full.model[[m]],new0.nm,type=type))/margin)
      else
        denm2<-cbind(denm2,(predict(full.model[[m]],new1.nm,best.iter1[m])-predict(full.model[[m]],new0.nm,best.iter1[m]))/margin)
      }
    denm1<-rbind(denm1,denm2)
    denm1.2=rbind(denm1.2,as.matrix(denm3))
    }
    denm1.0[[l]]<-denm1 
    denm1.1[[l]]<-denm1.2 
    te1.0[[l]]<-te1
  } 
  
  
  #4.4 get the indirect effects
  denm<-NULL
  denm1<-NULL
  te<-NULL
  ie<-NULL
  for (l in 1:nmod)
  {level=moder.level1$levels[,l]
   n.new=sum(level)
   denm[[l]]<-apply(denm1.0[[l]],2,col_mean,n.new)
   denm1[[l]]<-apply(denm1.1[[l]],2,col_mean,n.new)
   te0<-matrix(apply(te1.0[[l]],1,mean),n.new)
   colnames(te0)=paste("y",1:ncol(y),sep="")
   #te[[l]]<-te0
   temp1<-ncol(denm[[l]])/ncol(te0)
   temp2<-NULL
  for(temp in 1:temp1)
    temp2<-cbind(temp2,te0)
  ie[[l]]<-temp2-denm[[l]]
  ie[[l]][,1:ncol(y)]=matrix(rep(te0,ncol(y)),ncol=ncol(y))-denm1[[l]]      #the total indirect effect
  te[[l]]=as.matrix(ie[[l]][,1:ncol(y)]+denm[[l]][,1:ncol(y)])                    #the total effect
  colnames(te[[l]])=paste("y",1:ncol(y),sep="")
  if(!is.null(listm$multi)) 
    colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",colnames(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
  else 
    colnames(denm[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("de",colnames(x)[listm$single]),each=ncol(y)),sep=".")
  if(!is.null(listm$multi))
    colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",colnames(x)[listm$single],paste("j",1:listm$multi[[1]],sep="")),each=ncol(y)),sep=".")
  else 
    colnames(ie[[l]])<-paste(paste("y",1:ncol(y),sep=""),rep(c("all",colnames(x)[listm$single]),each=ncol(y)),sep=".")
  }
  names(te)<-paste(vari,moder.level1$moder.level,sep=".")
  names(denm)<-paste(vari,moder.level1$moder.level,sep=".")
  names(ie)<-paste(vari,moder.level1$moder.level,sep=".")
  
  a<-list(denm=denm,ie=ie,te=te,model=list(MART=nonlinear,Survival=surv, type=type, 
          model=full.model,best.iter=best.iter1),pred.new=pred.new,w.new=w.new,
          data=data,distmgivenx=distmgivenx,mod=T)
  class(a)<-"med"
  return(a)
}

anymissing<-function(vec)
{if(sum(is.na(vec))>0)
  return(F)
  else return(T)
}

mod.level<-function(vari=NULL,x=NULL,cova=NULL,continuous.resolution=10,w)
{pre=F
post=F
moder.level=NULL
moder=NULL
temp.q=NULL
if(is.null(w))
  w=rep(1,nrow(x))

if(sum(grep(vari,colnames(x)))>0)
{post=T #as a post moderator
moder=x[,vari]}
else if(sum(grep(vari,names(cova)))>0)
{pre=T  #as a pre moderator
moder=cova[,vari]}

temp.all1=NULL
if(is.factor(moder)){
  if(continuous.resolution==10)
    moder.level=levels(x[,vari])
  else
    moder.level=continuous.resolution
  for (le in moder.level)
    temp.all1=cbind(temp.all1,moder==le)}
else{
  if(length(continuous.resolution==1))
    temp.q=quantile(moder,probs=(seq(0,1,by=1/continuous.resolution))[-1],na.rm=T)
  else
    temp.q=continuous.resolution
  for(i in 1:length(temp.q))
  {if (i==1)
    temp.all=(moder<=temp.q[i])
  else
    temp.all=(moder<=temp.q[i] & moder>temp.q[i-1])
  temp.all[is.na(temp.all)]=F
  temp.all1=cbind(temp.all1,temp.all)
 # browser()
  moder.level=c(moder.level,weighted.mean(moder[temp.all],w[temp.all],na.rm=T))
  }
}
list(n.moder.level=length(moder.level),moder.level=moder.level,cont.moder.q=temp.q,moder=moder,levels=temp.all1)
}

data=mma1$data
x=data$x
y=data$y
dirx=data$dirx
contm=data$contm
catm=data$catm
binm=data$binm
jointm=data$jointm
cova=data$cova
allm=c(contm,catm)
xnames<-colnames(x)
pred_names<-colnames(dirx)
cova_names<-colnames(cova)
ynames=colnames(y)

surv=mma1$model$Survival
type=mma1$model$type
nonlinear=mma1$model$MART

if(is.null(x.new))
{x.new=x
pred.new=dirx
cova.new=cova
w.new=w}

if(!is.null(w.new)){
  if(is.null(cova.new))
    nonmissing1<-apply(cbind(pred.new,w.new),1,anymissing)
  else
    nonmissing1<-apply(cbind(cova.new,pred.new,w.new),1,anymissing)}
else{
  if(is.null(cova.new))
    nonmissing1<-apply(pred.new,1,anymissing)
  else
    nonmissing1<-apply(cbind(cova.new,pred.new),1,anymissing)}
x.new<-x.new[nonmissing1,]
colnames(x.new)=xnames
w.new<-w.new[nonmissing1]
pred.new<-data.frame(pred.new[nonmissing1,])
colnames(pred.new)<-pred_names
if(!is.null(cova.new))
{cova.new<-data.frame(cova.new[nonmissing1,])
colnames(cova.new)<-cova_names}

mod.level1=mod.level(vari,x.new,cova.new,continuous.resolution,w.new)

if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")


# if(ncol(x.new)>length(unique(contm,binm,catm)))
#  covay.new=x.new[,-unique(c(contm,binm,catm))]
# else covay.new=NULL

if(is.null(catm))
{multi=jointm
name1<-NULL                       #added in the new program
if (!is.null(multi))              #added in the new program, in case that multi is NULL
  name1<-paste("j",1:multi[[1]],sep="")}
else if(is.null(jointm))
{multi=catm
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,colnames(x)[multi[[i]][1]])}
else {temp1<-catm
temp2<-jointm
temp1[[1]]=catm[[1]]+jointm[[1]]
temp2[[1]]<-NULL
multi=append(temp1,temp2)
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,colnames(x)[multi[[i]][1]])
name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
listm=list(single=c(contm,binm),multi=multi)

ny=ncol(y)
nx=1#ncol(dirx)
nmod=mod.level1$n.moder.level
if(!is.null(mma1$all_boot))
  n2=nrow(mma1$all_boot)
te<-matrix(0,n2+1,ny*nmod)
de<-matrix(0,n2+1,ny*nmod)
mul<-ifelse(is.null(multi),0,multi[[1]])        #added in the new program, in case multi is null
ie<-matrix(0,n2,ny*(1+length(listm$single)+mul))   #added in the new program
ie1<-matrix(0,nmod,ny*(1+length(listm$single)+mul))   #added in the new program
if(!is.null(listm$multi))
{dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[listm$single],name1),each=ny),sep=".")
colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[listm$single],name1),each=ny),sep=".")
rownames(ie1)<-mod.level1$moder.level}
else 
{dimnames(ie)[[2]]<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[listm$single]),each=ny),sep=".")
colnames(ie1)<-paste(paste("y",1:ny,sep=""),rep(c("all",colnames(x)[listm$single]),each=ny),sep=".")
rownames(ie1)<-mod.level1$moder.level}
ie<-rep(list(ie),nmod)
names(ie)<-mod.level1$moder.level
#browser()

temp=mod.contx(vari,continuous.resolution,x,y,dirx,binm,contm,catm,jointm,cova, n=n,x.new=x.new,
               pred.new=pred.new, cova.new=cova.new, w=w, w.new=w.new,
               full.model=mma1$model$model,best.iter1=mma1$model$best.iter,surv=surv,type=type,
               moder.level1=mod.level1,nonlinear=nonlinear,df1=df1,n2,xj,xmod,margin=margin)
#temp=temp.med
for (l in 1:nmod)
  if(is.null(w.new))
  {ie1[l,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
   te[1,((l-1)*ny+1):(l*ny)]=apply(temp$te[[l]],2,mean,na.rm=T)
   de[1,((l-1)*ny+1):(l*ny)]=apply(as.matrix(temp$denm[[l]][,((l-1)*ny+1):(l*ny)]),2,mean,na.rm=T)
  }
  else
  {level=mod.level1$levels[,l]
   te[1,((l-1)*ny+1):(l*ny)]<-apply(temp$te[[l]],2,weighted.mean,na.rm=T,w=w.new[level])
   de[1,((l-1)*ny+1):(l*ny)]<-apply(as.matrix(temp$denm[[l]][,1:ny]),2,weighted.mean,na.rm=T,w=w.new[level]) 
   ie1[l,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new[level])  #first row is the estimated value
  }

te1<-rep(list(NULL),nmod)                      #to store the mediation effects on predictor
de1<-rep(list(NULL),nmod)
ie2<-rep(list(NULL),nmod)
names(ie2)<-mod.level1$moder.level
names(te1)<-mod.level1$moder.level
names(de1)<-mod.level1$moder.level

if(!is.null(mma1$all_boot)){
  n2=nrow(mma1$all_boot)
  for (i in 1:n2)
  {boots<-mma1$all_boot[i,]
  x1<-data.frame(x[boots,])
  colnames(x1)=xnames
  y1<-data.frame(y[boots,])
  colnames(y1)=ynames
  dirx1<-data.frame(dirx[boots,])
  colnames(dirx1)=pred_names
  if(!is.null(cova))
    {if(length(grep("for.m",names(cova)))==0)
  {cova1<-data.frame(cova[boots,])
  names(cova1)=cova_names}
  else if(!is.null(cova))
  {cova1=cova
  cova1[[1]]=data.frame(cova[[1]][boots,])
  names(cova1[[1]])=cova_names}}
  else
    cova1=NULL
  
  
  temp<-mod.contx(vari=vari,continuous.resolution=continuous.resolution,x=x1,y=y1,
                  dirx=dirx1,binm=binm,contm=contm,catm=catm,jointm=jointm,cova=cova1, n=n,x.new=x.new,
                  pred.new=pred.new, cova.new=cova.new, w=NULL, w.new=w.new,
                  full.model=mma1$all_model[[i]],best.iter1=mma1$all_iter[i,],surv=surv,type=type,
                  moder.level1=mod.level1,nonlinear=nonlinear,df1=df1,n2,xj,xmod,margin=margin)
  
  for (l in 1:nmod)
  {if(is.null(w.new))
   {te[1+i,((l-1)*ny+1):(l*ny)]<-apply(temp$te[[l]],2,mean,na.rm=T)
    de[1+i,((l-1)*ny+1):(l*ny)]<-apply(as.matrix(temp$denm[[l]][,((l-1)*ny+1):(l*ny)]),2,mean,na.rm=T)
    ie[[l]][i,]<-apply(temp$ie[[l]],2,mean,na.rm=T)  #first row is the estimated value
   }
   else
   {level=mod.level1$levels[,l]
    te[1+i,((l-1)*ny+1):(l*ny)]<-apply(temp$te[[l]],2,weighted.mean,na.rm=T,w=w.new[level])
    de[1+i,((l-1)*ny+1):(l*ny)]<-apply(as.matrix(temp$denm[[l]][,1:ny]),2,weighted.mean,na.rm=T,w=w.new[level])
    ie[[l]][i,]<-apply(temp$ie[[l]],2,weighted.mean,na.rm=T,w=w.new[level])  #first row is the estimated value
   }
  
   te1[[l]]<-cbind(te1[[l]],temp$te[[l]])
   de1[[l]]<-cbind(de1[[l]],as.matrix(temp$denm[[l]][,1:ny]))
   ie2[[l]]<-rbind(ie2[[l]],temp$ie[[l]])
  }
  
  print(i)
  }
  
  colnames(te)<-paste(paste("y",1:ncol(y),sep=""),rep(mod.level1$moder.level,each=ncol(y)),sep=".")
  colnames(de)<-paste(paste("y",1:ncol(y),sep=""),rep(mod.level1$moder.level,each=ncol(y)),sep=".")
}

a<-list(estimation=list(ie=ie1,te=te[1,],de=de[1,]),bootsresults=list(ie=ie,te=te[-1,],de=de[-1,]),model=mma1$model,
        data=list(x=x,y=y,dirx=dirx,binm=binm,contm=contm,catm=catm, jointm=jointm, cova=cova, binpred=F),
        boot.detail=list(pred.new=pred.new,cova.new=cova.new,te1=te1,de1=de1,ie1=ie2),w.new=w.new, pred.new=pred.new,
        moder.level=mod.level1,mod=T)
class(a)<-"mma"
return(a)
}

if(!is.null(mma1$a.binx))
{binpred=mma1$a.binx$data$binpred
 contpred=mma1$a.binx$data$contpred 
 catpred=mma1$a.binx$data$catpred 
}
else
{binpred=mma1$a.contx$data$binpred
 contpred=mma1$a.contx$data$contpred 
 catpred=mma1$a.contx$data$catpred 
}

a.binx=NULL
a.contx=NULL

if(xj%in%contpred)
  {if(is.null(w))
    w=rep(1,nrow(mma1$a.contx$data$x))
   mma1$a.binx$data$binpred=F
   a.contx<-boot.mod.contx(mma1$a.contx,vari,continuous.resolution=continuous.resolution,
                    w=w,n=n,x.new=x.new,w.new=w.new,pred.new=pred.new,
                    cova.new=cova.new,xj=xj,df1=df1,xmod=xmod,margin=margin)
}
else if(xj%in%binpred)
{if(is.null(w))
  w=rep(1,nrow(mma1$a.binx$data$x))
 mma1$a.binx$data$binpred=T
 a.binx<-boot.mod.binx(mma1$a.binx,vari,continuous.resolution=continuous.resolution,n=n,w=w,xj=xj,xmod=xmod)
}
else
{z11=rep(F,length(catpred))
  for (i in 1:length(catpred))
   z11[i]=xj%in%catpred[[i]]
 i=(1:length(catpred))[z11]
 if(is.null(w))
  w=rep(1,nrow(mma1$a.binx$data$x))
 mma1$a.binx$data$binpred=T
 a.binx<-boot.mod.binx(mma1$a.binx,vari,continuous.resolution=continuous.resolution,n=n,w=w,xj=catpred[[i]],xmod=xmod)
}
  
a<-list(a.binx=a.binx,a.contx=a.contx,pred=list(binpred=binpred,catpred=catpred,contpred=contpred))
class(a)="mma"
return(a)
}


plot2.mma<-function(x,...,vari,xlim=NULL,alpha=0.95,quantile=F,moderator,xj=1)
{plot2.temp<-function(x,...,vari,xlim=NULL,alpha=0.95,quantile=F,moderator,xj=1){
  marg.den<-function(x,y,w=NULL) #added w
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

plot_ci<-function(df1,xlab="x",ylab="IE",sub=NULL)
{ plot(df1$x, df1$F, ylim = range(c(df1$L,df1$U),na.rm=T), type = "l",xlab=xlab,ylab=ylab,sub=sub)
  polygon(c(df1$x,rev(df1$x)),c(df1$L,rev(df1$U)),col = "grey75", border = FALSE)
  lines(df1$x, df1$F, lwd = 2)
  lines(df1$x, df1$U, col="red",lty=2)
  lines(df1$x, df1$L, col="red",lty=2)}

nx<-ncol(x$data$dirx)
ny<-ncol(x$data$y)
nmod=length(x$moder.level$moder.level)
op <- par(no.readonly = TRUE) # the whole list of settable par's.
data=x$data
mname<-ifelse(is.character(vari),vari,names(data$x)[vari])
vari=mname
if(is.null(xlim) & !is.factor(x$data$x[,grep(vari,names(x$data$x))]))
  xlim=range(x$data$x[,grep(vari,colnames(x$data$x))],na.rm=T)

if (x$model[1]==T) 
  for (m in 1:ny) {
    full.model=x$model$model[[m]]
    best.iter=x$model$best.iter[m]
    if(data$binpred)
    {d<-rep(0,nrow(data$dirx))
    for(l in 1:nx)
      d[data$dirx[,l]==1]<-l
    if(!is.factor(data$x[,vari]))
    { par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      if(full.model$distribution=="gaussian")
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,xlim=xlim)))
      else if(full.model$distribution=="coxph")
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator),xlim=xlim)))
      else
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,xlim=xlim,type="response")))

      par(mfrow=c(max(2,min(5,ceiling(nmod/2))),nx+1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      for(q1 in 1:nmod){
        if(is.factor(x$moder.level$moder))
          temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
        else if (q1==1)
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & !is.na(x$moder.level$moder))
        else
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & x$moder.level$moder>x$moder.level$cont.moder.q[q1-1] & !is.na(x$moder.level$moder))
        
        overlapHist(a=data$x[temp.all,vari],b=as.matrix(d[temp.all]),xlim=xlim,xname=paste(moderator, "=", x$moder.level$moder.level[q1],", Predictor"),w=data$w[temp.all])} # added w
    }
    else{
      if(full.model$distribution=="gaussian")
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator),best.iter)))
      else if(full.model$distribution=="coxph")
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator))))
      else
        suppressWarnings(print(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,type="response")))
      par(mfrow=c(max(2,min(5,ceiling(nmod/2))),nx+1),mar=c(5,5,1,1),oma=c(3,2,5,4)) 
      for(q1 in 1:nmod){
        if(is.factor(x$moder.level$moder))
          temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
        else if (q1==1)
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & !is.na(x$moder.level$moder))
        else
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & x$moder.level$moder>x$moder.level$cont.moder.q[q1-1] & !is.na(x$moder.level$moder))
       
        temp1<-NULL
        if (is.null(data$w)) #
        {temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0 & temp.all,vari])))
        for (j in 1:nx)
          temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1 & temp.all,vari])))
        barplot(prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0 & temp.all,vari])),ylim=c(0,max(temp1,na.rm=T)),
                ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1], ", Predictor at the Reference Level: pred=",0,sep=""))     
        #browser()
        for (j in 1:nx)
          barplot(prop.table(table(data$x[data$dirx[,j]==1 & temp.all,vari])),ylim=c(0,max(temp1,na.rm=T)),
                  ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1], colnames(data$dirx)[j], ", pred=",j,sep=""))}
        else #
        {temp1<-c(temp1,weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0 & temp.all,vari],data$w))#
        for (j in 1:nx) #
          temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1 & temp.all,vari],data$w))#
        barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0 & temp.all,vari]),ylim=c(0,max(temp1,na.rm=T)),#
                ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1], ", Predictor at the Reference Level, pred=", j,sep=""))
        for (j in 1:nx)#
          barplot(weighted.prop.table(data$x[data$dirx[,j]==1 & temp.all,vari]),ylim=c(0,max(temp1,na.rm=T)),#
                  ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1], colnames(data$dirx)[j], ", pred=", j,sep=""))} #
      }}
    }
    else
    {par(mfrow=c(ceiling(nmod/2),2),mar=c(5,5,1,1),oma=c(3,2,5,4))
      for (l in 1:nmod)
      {temp.ie.detail<-as.matrix(x$boot.detail$ie1[[l]][,grep(mname,colnames(x$boot.detail$ie1[[l]]))])  #
       level=x$moder.level$levels[,l]
      ie1<-boot.ci(x$pred.new[level,xj],matrix(temp.ie.detail[,m],nrow=sum(level)),alpha,quantile)
      plot_ci(ie1,xlab=names(data$dirx)[xj],ylab=paste("IE on",colnames(data$y)[m]),sub=x$moder.level$moder.level[l])}
      par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))   
      if(!is.factor(data$x[,vari]))
      {if(full.model$distribution=="gaussian")
        print(suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,xlim=xlim)))
        else if(full.model$distribution=="coxph")
          print(suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator),xlim=xlim)))
        else
          print(suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,xlim=xlim,type="response")))
#        if(nx>1)
#          for (i in 1:(nx-1))
#            plot(1, type="n", axes=F, xlab="", ylab="")
        par(mfrow=c(ceiling(nmod/2),2),mar=c(5,5,1,1),oma=c(3,2,5,4))
        for(l in 1:nmod){
       #browser()
          if(is.factor(x$moder.level$moder))
            temp.all=(x$moder.level$moder==x$moder.level$moder.level[l] & !is.na(x$moder.level$moder))
          else
            temp.all=x$moder.level$levels[,l]
          
          a<-marg.den(data$dirx[temp.all,xj],data$x[temp.all,vari],data$w[temp.all]) #added data$w
          scatter.smooth(a[,1],a[,2],family="gaussian",xlab=colnames(data$dirx)[xj],ylim=xlim,ylab=paste("Mean",mname,sep="."),sub=paste(moderator, "at", x$moder.level$moder.level[l]))
          axis(1,at=data$x[temp.all,vari],labels=F)}
      }
      else
      {par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
        if(full.model$distribution=="gaussian")
          suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator),best.iter))
        else if(full.model$distribution=="coxph")
          suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator)))
        else
          suppressWarnings(plot.gbm(full.model, i.var=c(vari,moderator),best.iter,type="response"))
        par(mfrow=c(ceiling(nmod/2),2),mar=c(5,5,1,1),oma=c(3,2,5,4))
        #        if(nx>1)
        #          for (i in 1:(nx-1))
        #            plot(1, type="n", axes=F, xlab="", ylab="")
        for(l in 1:nmod){
          if(is.factor(x$moder.level$moder))
            temp.all=(x$moder.level$moder==x$moder.level$moder.level[l] & !is.na(x$moder.level$moder))
          else if (l==1)
            temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & !is.na(x$moder.level$moder))
          else
            temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & x$moder.level$moder>x$moder.level$cont.moder.q[l-1] & !is.na(x$moder.level$moder))
          
          plot(data$x[temp.all,vari],data$dirx[temp.all,xj],ylab=colnames(data$dirx)[l],xlab="",sub=paste(moderator,"at", x$moder.level$moder.level[l]))}}
    }
  }
else
  for (m in 1:ny) 
  {full.model=x$model$model[[m]]
  coef<-full.model$coefficients[grep(vari,names(full.model$coefficients))] #plot the straight line instead of the loess line
  if(is.null(full.model$na.action))
    data1<-full.model$data[,vari]
  else
  {data1<-full.model$data[-full.model$na.action,vari]
  x$moder.level$moder=x$moder.level$moder[-full.model$na.action]}
  if(x$model$Survival[m] & is.null(x$model$best.iter)) #for cox model
    data1<-x$data$x[,vari]
  
  if(data$binpred)
  {d<-rep(0,nrow(data$dirx))
  for(l in 1:nx)
    d[data$dirx[,l]==1]<-l
  if(!is.factor(data$x[,vari]))
  {par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    if(!x$model$Survival[m])
      b1<-full.model$family$linkfun(full.model$fitted.values) #added data$w
    else
      b1<-predict(full.model,type=x$model$type)  #added data$w
    plot(data1,b1,type="n",xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
    for(l in 1:nmod){
      if(is.factor(x$moder.level$moder))
        temp.all=(x$moder.level$moder==x$moder.level$moder.level[l] & !is.na(x$moder.level$moder))
      else if (l==1)
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & !is.na(x$moder.level$moder))
      else
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & x$moder.level$moder>x$moder.level$cont.moder.q[l-1] & !is.na(x$moder.level$moder))
      
      b<-marg.den(data1[temp.all],b1[temp.all],data$w[temp.all])  #added data$w
      #browser()      
      points(b,col=l)
      if(length(coef)>1)
      {b2=coef[grep(x$moder.level$moder.level[l],names(coef))]+coef[vari]
      b3=ifelse(is.factor(x$moder.level$moder.level),1,x$moder.level$moder.level[l])
      b2.1=coef[grep(x$moder.level$moder.level[l],names(coef))]*b3+coef[vari]*mean(b[,1])
      abline(a=mean(b[,2],na.rm=T)-b2.1,b=b2,col=l)}
      else
        abline(a=mean(b[,2],na.rm=T)-coef[vari]*mean(b[,1]),b=coef[vari],col=l)
      
      axis(1,at=data1,labels=F)}
    par(mfrow=c(max(2,min(5,ceiling(nmod/2))),nx+1),mar=c(5,5,1,1),oma=c(3,2,5,4))  
    for(q1 in 1:nmod){
      if(is.factor(x$moder.level$moder))
        temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
      else if (q1==1)
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & !is.na(x$moder.level$moder))
      else
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & x$moder.level$moder>x$moder.level$cont.moder.q[q1-1] & !is.na(x$moder.level$moder))
      overlapHist(a=data$x[temp.all,vari],b=as.matrix(d[temp.all]),xlim=xlim,xname=paste(moderator, "=", x$moder.level$moder.level[q1],", Predictor"),data$w[temp.all])}  #added data$w
  }
  else{
    if(is.factor(x$moder.level$moder))
    {if(!x$model$Survival[m])
      b1=full.model$fitted.values
    else
      b1=predict(full.model,se.fit=T,type=x$model$type)$fit
    print(xyplot(b1~data1|x$moder.level$moder,ylab=paste("f(",mname,")",sep=""),xlab=mname))
    }
    else
    {if(!x$model$Survival[m])
      print(levelplot(full.model$fitted.values~data1*x$moder.level$moder,ylab=moderator,xlab=mname))
      else
        print(levelplot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1*x$moder.level$mode,ylab=moderator,xlab=mname))}
    
    temp1<-NULL
    if(is.null(data$w)){ #
      temp1<-c(temp1,prop.table(table(data$x[apply(data$dirx!=0,1,sum)==0,vari])))
      for (j in 1:ncol(data$dirx))
        temp1<-c(temp1,prop.table(table(data$x[data$dirx[,j]==1,vari])))
      par(mfrow=c(max(2,min(5,ceiling(nmod/2))),nx+1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      for(q1 in 1:nmod){
        if(is.factor(x$moder.level$moder))
          temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
        else if (q1==1)
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & !is.na(x$moder.level$moder))
        else
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & x$moder.level$moder>x$moder.level$cont.moder.q[q1-1] & !is.na(x$moder.level$moder))
        
        barplot(prop.table(table(data$x[apply(data$dirx!=0 & temp.all,1,sum)==0,vari])),ylim=c(0,max(temp1,na.rm=T)),
                ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1],"Predictor at the reference level"))
        for (j in 1:ncol(data$dirx))
          barplot(prop.table(table(data$x[data$dirx[,j]==1 & temp.all,vari])),ylim=c(0,max(temp1,na.rm=T)),
                  ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1],"Predictor at",colnames(data$dirx)[j]))}
    }#
    else#
    {temp1<-c(temp1,weighted.prop.table(table(data$x[apply(data$dirx,1,sum)==0,grep(vari,names(data$x))],data$w)))
    for (j in 1:ncol(data$dirx))#
      temp1<-c(temp1,weighted.prop.table(data$x[data$dirx[,j]==1,grep(vari,names(data$x))],data$w))#
    par(mfrow=c(max(2,min(5,ceiling(nmod/2))),nx+1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    for(q1 in 1:nmod){
      if(is.factor(x$moder.level$moder))
        temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
      else if (q1==1)
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & !is.na(x$moder.level$moder))
      else
        temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[q1] & x$moder.level$moder>x$moder.level$cont.moder.q[q1-1] & !is.na(x$moder.level$moder))
      barplot(weighted.prop.table(data$x[apply(data$dirx!=0,1,sum)==0 & temp.all,vari],data$w[temp.all]),ylim=c(0,max(temp1)),#
              ylab="Prop",sub=paste(moderator, "=", x$moder.level$moder.level[q1],"Predictor at the reference level")) #
      for (j in 1:ncol(data$dirx))#
        barplot(weighted.prop.table(data$x[data$dirx[,j]==1 & temp.all,vari],data$w[temp.all]),ylim=c(0,max(temp1)),#
                ylab="Prop",sub=colnames(data$dirx)[j])} #
    }}
  }
  else
  {par(mfrow=c(ceiling(nmod/2),2),mar=c(5,5,1,1),oma=c(3,2,5,4))
    for (l in 1:nmod) {
      level=x$moder.level$levels[,l]
      temp.ie.detail<-as.matrix(x$boot.detail$ie1[[l]][,grep(mname,colnames(x$boot.detail$ie1[[l]]))])  #
      ie1<-boot.ci(x$pred.new[level,xj],matrix(temp.ie.detail[,m],nrow=sum(level)),alpha,quantile)
      plot_ci(ie1,xlab=colnames(data$dirx)[xj],sub=x$moder.level$moder.level[l])}
    
    if(!is.factor(data$x[,vari]))
    {if(!x$model$Survival[m])
      b1<-full.model$family$linkfun(full.model$fitted.values) #added data$w
    else
      b1<-predict(full.model,se.fit=T,type=x$model$type)$fit #added data$w
    
    #par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
    
    par(mfrow=c(max(2,min(5,ceiling(nmod/2))),2),mar=c(5,5,1,1),oma=c(3,2,5,4))
    for(q1 in 1:nmod){ 
      if(is.factor(x$moder.level$moder))
        temp.all=(x$moder.level$moder==x$moder.level$moder.level[q1] & !is.na(x$moder.level$moder))
      else 
        temp.all=x$moder.level$levels[,q1]
      b<-marg.den(data1[temp.all],b1[temp.all],data$w[temp.all]) #added data$w
      plot(data1,b1,type="n",xlab=mname,ylab=paste("f(",mname,")",sep=""),xlim=xlim)
      
      points(b,col=q1)
      if(length(coef)>1){#browser()
        b3=ifelse(is.factor(x$moder.level$moder.level),1,x$moder.level$moder.level[q1])
        if(is.factor(x$moder.level$moder))
         b2=coef[grep(x$moder.level$moder.level[q1],names(coef))]+coef[vari]
        else
          b2=coef[intersect(grep(vari,names(coef)),grep(moderator,names(coef)))]*b3
        
        if(is.factor(x$moder.level$moder))
          b2.1=coef[grep(x$moder.level$moder.level[q1],names(coef))]*b3+coef[vari]*mean(b[,1])
        else
          b2.1=coef[intersect(grep(vari,names(coef)),grep(moderator,names(coef)))]*b3+coef[vari]*mean(b[,1])
        abline(a=mean(b[,2],na.rm=T)-b2.1,b=b2,col=q1)} #
      else
        abline(a=mean(b[,2],na.rm=T)-coef*mean(b[,1]),b=coef,col=q1)
      # browser() 
      axis(1,at=data1,labels=F)
      #if(nx>1)
      #  for (i in 1:(nx-1))
      #    plot(1, type="n", axes=F, xlab="", ylab="")
      #for(l in 1:nx){
      a<-marg.den(x$pred.new[temp.all,xj],data$x[temp.all,vari],data$w[temp.all])   #added data$w
      scatter.smooth(a[,1],a[,2],family="gaussian", xlab=colnames(data$dirx)[xj],ylim=xlim,ylab=paste("Mean",mname,sep="."),
                     sub=x$moder.level$moder.level[q1])}
    }
    else
    {par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(3,2,5,4))
      if(is.factor(x$moder.level$moder))
      {if (!x$model$Survival[m])
        print(xyplot(full.model$fitted.values~data1|x$moder.level$moder,ylab=paste("f(",mname,")",sep=""),xlab=mname))
        else
          print(xyplot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1|x$moder.level$mode,ylab=paste("f(",mname,")",sep=""),xlab=mname))}
      else
      {if (!x$model$Survival[m])
        print(levelplot(full.model$fitted.values~data1*x$moder.level$moder,ylab=moderator,xlab=mname))
        else
          print(levelplot(predict(full.model,se.fit=T,type=x$model$type)$fit~data1*x$moder.level$mode,ylab=moderator,xlab=mname))}
      # if(nx>1)
      #    for (i in 1:(nx-1))
      #     plot(1, type="n", axes=F, xlab="", ylab="")
      for(l in 1:nmod){
        if(is.factor(x$moder.level$moder))
          temp.all=(x$moder.level$moder==x$moder.level$moder.level[l] & !is.na(x$moder.level$moder))
        else if (l==1)
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & !is.na(x$moder.level$moder))
        else
          temp.all=(x$moder.level$moder<=x$moder.level$cont.moder.q[l] & x$moder.level$moder>x$moder.level$cont.moder.q[l-1] & !is.na(x$moder.level$moder))
        plot(data$x[temp.all,vari],x$moder.level$pred.uniq[[l]],ylab=colnames(data$dirx)[xj],xlab="",sub=paste(moderate, "at",x$moder.level$moder.level[l]))}}
  }
  }
par(op)
}
contpred=x$pred$contpred
catpred=x$pred$catpred
binpred=x$pred$binpred


if(xj%in%contpred)
  plot2.temp(x=x$a.contx,vari=vari,xlim=xlim,alpha=alpha,quantile=quantile,moderator=moderator,xj=xj)
else if(xj%in%binpred)
  plot2.temp(x=x$a.binx,vari=vari,xlim=xlim,alpha=alpha,quantile=quantile,moderator=moderator,xj=xj)
else
{z11=rep(F,length(catpred))
 for (i in 1:length(catpred))
  z11[i]=xj%in%catpred[[i]]
 i=(1:length(catpred))[z11]
 plot2.temp(x=x$a.binx,vari=vari,xlim=xlim,alpha=alpha,quantile=quantile,moderator=moderator,xj=catpred[[i]])
}
}



