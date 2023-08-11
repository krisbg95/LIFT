library(DataExplorer) 
library(readr) 
library(corrplot) 
library(Hmisc)

data=read.csv("BRAF.csv",row.names=1)  
data[,1] = NULL 
colnames (data) = c("Physical","Living", "Cognitive", "Emotional")

# mcorr <- rcorr(as.matrix(data),type="spearman")
# mcol = colorRampPalette(c("blue","red"))(200) 
# corrplot(mcorr$r, method="number",col=mcol,type="lower", order="hclust",p.mat = mcorr$P, sig.level = 0.05, insig = "blank",tl.col = "black")

pairs.panels(data, 
             method = "spearman",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = FALSE,lm = T,stars = T,ci = T) 

dev.print(svg,"BRAF_Comparisons.svg",height=8,width=8)
