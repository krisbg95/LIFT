#------------------------------------------------------------------------------------------------------
# install.packages("mmabig")
library(mmabig)
setwd("./")

medvars <- read.table("./varsSC.txt", quote="\"", comment.char="")
pred <- read.table("./varspred.txt", quote="\"", comment.char="",)
pred <- data.frame(pred[-19,])
y <- read.table("./varsy.txt", quote="\"", comment.char="")
y <- data.frame(y[-c(19),])
cova <- read.table("./varscova.txt", quote="\"", comment.char="")
cova <- cova[-c(19),]
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

g32_ind<-pred$Group!=1
g32_cova<-cova[g32_ind,]
g32_medvars<-medvars[g32_ind,]
g32_pred<-pred[g32_ind,]
g32_y<-y[g32_ind,]

cova=data.frame(g32_cova)
mediator=c(1:ncol(g32_medvars))

x<-data.frame(cbind(g32_medvars,g32_cova))

data.e1<-data.org.big(x,y=data.frame(g32_y),mediator=mediator,pred=data.frame(g32_pred),predref="2",alpha=1,alpha1=0.1,alpha2=0.1,testtype=2)
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
frmsig3<-as.data.frame(frm[,c])
med_list3<-NULL

for(r in 1:(ncol(frm)/2)){
  val1g3<-frm[4,r+(ncol(frm)/2)]
  val2g3<-frm[5,r+(ncol(frm)/2)]
  mulvalg3<-val1g3*val2g3
  if (mulvalg3>0){
    frmsig3[,c]<-frm[,r]
    tempnames<-colnames(frm)[r]
    med_list3<-cbind(med_list3,tempnames)
    c=c+1}}

colnames(frmsig3)<-med_list3
row.names(frmsig3)<-row.names(frm)

write.table(data.frame(colnames(frmsig3)), 'med_list3SC_g32.csv', append=T, sep=',')
write.table(data.frame(frmsig3), 'SC3_g32.csv', append=T, sep=',')

save.image("./workspace_SC_g32.RData")
