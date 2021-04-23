#------------------------------------------------------------------------------------------------------
# install.packages("mmabig")
library(mmabig)
setwd("./")

medvars <- read.table("./varsrsFC.txt", quote="\"", comment.char="")
pred <- read.table("./varspred.txt", quote="\"", comment.char="")
y <- read.table("./varsy.txt", quote="\"", comment.char="")
cova <- read.table("./varscova.txt", quote="\"", comment.char="")
medheader <- read.table("./headervars.txt", sep="\t", quote="\"", comment.char="")
medheader <- medheader[1:(length(medheader)-1)]
covaheader <- read.table("./headercova.txt", sep="\t", quote="\"", comment.char="")
covaheader <- covaheader[1:(length(covaheader)-1)]
# colnames(medvars) <- medheader
colnames(cova) <- covaheader
colnames(pred) = 'Group'
colnames(y) = 'Chalderdiff'

pred$Group <- as.factor(pred$Group)
cova$CENTRE <- as.factor(cova$CENTRE)
cova$SEX <- as.factor(cova$SEX)

g12_ind<-pred$Group!=3
g12_cova<-cova[g12_ind,]
g12_medvars<-medvars[g12_ind,]
g12_pred<-pred[g12_ind,]
g12_y<-y[g12_ind,]

cova=data.frame(g12_cova)
mediator=c(1:ncol(g12_medvars))

x<-data.frame(cbind(g12_medvars,g12_cova))

data.e1<-data.org.big(x,y=data.frame(g12_y),mediator=mediator,pred=data.frame(g12_pred),predref="2",alpha=1,alpha1=0.1,alpha2=0.1,testtype=2)
# summary(data.e1,only=TRUE)

med.e1<-med.big(data.e1)
# med.e1

mma.e1<-mma.big(data=data.e1,alpha=1,alpha1=0.1,alpha2=0.1,n2=1000) # use only the test results.

alpha<-0.05
a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
x<-mma.e1
ny<-ncol(x$data$y)
nx<-ncol(x$data$dirx)
temp1<-x$bootresults
temp.0<-print(x$results,print=F)
temp2<-temp1
for(l in 1:ny)
  temp2[[l]]=temp2[[l]]%*%diag(1/temp2[[l]][nrow(temp2[[l]]),])
temp<-vector("list",ny)
names(temp)<-colnames(x$data$y)
temp.4<-vector("list",ny)                #save relative effects
names(temp.4)<-colnames(x$data$y)

for(l in 1:ny){
  temp.1<-NULL
  for (j in 1:nrow(temp1[[l]]))
    temp.1<-rbind(temp.1,matrix(temp1[[l]][j,],nrow=nx))
  temp[[l]]<-rbind(est=as.vector(temp.0$results[[l]]),
                   mean=as.vector(t(matrix(apply(temp.1,1,mean,na.rm=T),nrow=nx))),
                   sd=as.vector(t(matrix(apply(temp.1,1,sd,na.rm=T),nrow=nx))))
  temp[[l]]=rbind(temp[[l]],CI_upbd=temp[[l]][2,]+b2*temp[[l]][3,]/sqrt(1000))
  temp[[l]]=rbind(temp[[l]],CI_lwbd=temp[[l]][2,]+b1*temp[[l]][3,]/sqrt(1000))
  colnames(temp[[l]])=paste(rep(colnames(x$data$dirx),each=nrow(temp1[[l]])),rownames(temp1[[l]]),sep=".")
}

frm<-as.data.frame(temp)

c=1
frmsig1<-as.data.frame(frm[,c])
med_list1<-NULL

for(r in 1:(ncol(frm)/2)){
  val1g1<-frm[4,r]
  val2g1<-frm[5,r]
  mulvalg1<-val1g1*val2g1
  if (mulvalg1>0){
    frmsig1[,c]<-frm[,r]
    tempnames<-colnames(frm)[r]
    med_list1<-cbind(med_list1,tempnames)
    c=c+1}}

colnames(frmsig1)<-med_list1
row.names(frmsig1)<-row.names(frm)

write.table(data.frame(colnames(frmsig1)), 'med_list1rsFC_g12.csv', append=T, sep=',')
write.table(data.frame(frmsig1), 'rsFC1_g12.csv', append=T, sep=',')

save.image("./workspace_rsFC_g12.RData")
