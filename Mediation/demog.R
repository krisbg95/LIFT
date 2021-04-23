g1<-cova[pred$Group==1,]
g2<-cova[pred$Group==2,]
g3<-cova[pred$Group==3,]
mean(g1$AGE1)
sd(g1$AGE1)
mean(g2$AGE1)
sd(g2$AGE1)
mean(g3$AGE1)
sd(g3$AGE1)
count(g1,vars = 'SEX')
count(g1,vars = 'CENTRE')
count(g2,vars = 'SEX')
count(g2,vars = 'CENTRE')
count(g3,vars = 'SEX')
count(g3,vars = 'CENTRE')