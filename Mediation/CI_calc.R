
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


