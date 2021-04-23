#------------------------------------------------------------------------------------------------------
# install.packages("mmabig")
library(mmabig)

medvars <- read.table("./varsvol.txt", quote="\"", comment.char="")
pred <- read.table("./varspred.txt", quote="\"", comment.char="")
y <- read.table("./varsy.txt", quote="\"", comment.char="")
cova <- read.table("./varscova.txt", quote="\"", comment.char="")
medheader <- read.table("./headervars.txt", sep="\t", quote="\"", comment.char="")
medheader <- medheader[1:(length(medheader)-1)]
covaheader <- read.table("./headercova.txt", sep="\t", quote="\"", comment.char="")
covaheader <- covaheader[1:(length(covaheader)-1)]
colnames(medvars) <- medheader
colnames(cova) <- covaheader
colnames(pred) = 'Group'
colnames(y) = 'Chalderdiff'

pred$Group <- as.factor(pred$Group)
cova$CENTRE <- as.factor(cova$CENTRE)
cova$SEX <- as.factor(cova$SEX)

cova=data.frame(cova)
mediator=c(1:84)

x<-data.frame(cbind(medvars,cova))

# data.e1<-data.org.big(x,y=data.frame(y),mediator=mediator,pred=data.frame(pred),predref="2",testtype = 1)
# summary(data.e1,only=TRUE)
# 
# data.e1.2<-data.org.big(x,y=data.frame(y),mediator=mediator,pred=data.frame(pred),predref="2")
# summary(data.e1.2,only=TRUE)
# 
# med.e1<-med.big(data.e1)
# med.e1
# med.e1.2<-med.big(data.e1.2)
# med.e1.2
# 
# mma.e1<-mma.big()

data.e1<-data.org.big(x,y=data.frame(y),mediator=mediator,pred=data.frame(pred),predref="2",alpha=0.1,alpha2=0.1,testtype=1)
summary(data.e1,only=TRUE)

med.e1<-med.big(data.e1)
med.e1

mma.e1<-mma.big(data=data.e1,alpha=0.05,alpha2=0.05,n2=1000) # use only the test results.
summary(mma.e1)
# plot(mma.e1.2,vari="")

# medres<-med.big(data=medres)
# summary(medres)
# med_1.boot<-boot.med(data=med,n=100,n2=1000)
# 
# summary(med_1.boot,RE=F,alpha=0.05)
# plot(medres,RE=TRUE)
