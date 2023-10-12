## 0. Libraries
library(reshape2)  
library(effectsize)
library(ggplot2)
library(ggpubr)    
library(patchwork)
library(tidyverse) 
library(rcompanion) 
library(FSA) 
library(rstatix)

## 1. Upload data
load("All_scripts.RData") 
mclin = final$`Sequential/Clinical`$Fsel009
mfunctional = final$`Sequential/FC`$Fsel009 
mmulti = final$`Sequential/Multimodal`$Fsel010 
mSC = final$`Sequential/SC`$Fsel012 
mstruc = final$`Sequential/Structural`$Fsel007  

## 2. Structure data
listperf = list(mclin,mfunctional,mmulti,mSC,mstruc)  
names(listperf) = c("Clinical","FC","Multimodal","SC","Morphometric") 
listperfRMSE = lapply(listperf,function(x){x[,grepl("RMSE",colnames(x))]}) 
RMSE = lapply(listperfRMSE,function(x){x[,c(5,8)]}) 
RMSE = as.data.frame(unlist(RMSE,recursive = F))
RMSE$Morphometric.gaussprRadialRMSE = listperf$Morphometric$rfRMSE  
colnames(RMSE) = c("Clinical_gaussprRadial","Clinical_Base","FC_gaussprRadial","FC_Base","Multimodal_gaussprRadial","Multimodal_Base","SC_gaussprRadial","SC_Base","Morphometric_RF","Morphometric_Base")
# write.csv(RMSE,"Best_model_Base_RMSE.csv",row.names = F)
RMSE_long = melt(RMSE)
RMSEcompare = RMSE[,c(1,3,5,7,9)]
RMSEcompare_long = melt(RMSEcompare)
# write.csv(RMSEcompare_long,"Best_model_RMSE.csv",row.names = F)


## 3. Setup for baseline comparisons 
p_values <- numeric(5)
test_statistic <- numeric(5)
effect_size <- numeric(5)
n = nrow(RMSE)
for (i in 1:5) {
  wilcox_result <- wilcox.test(RMSE[, i*2-1], RMSE[, i*2], paired = TRUE)
  p_values[i] <- wilcox_result$p.value
  test_statistic[i] <- wilcox_result$statistic
  effect_size[i] <- effectsize::rank_biserial(RMSE[, i*2-1], RMSE[, i*2], paired = TRUE)[,1]
}


# p <- ggplot(RMSE_long, aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() + 
#   labs(x = "Model", y = "Performance") +
#   scale_fill_discrete(name = "variable") +
#   ggtitle("Model vs. Base Performance") +
#   theme_minimal()+ 
#   theme(legend.position = "none") 
# 
# com = levels(RMSE_long$variable) 
# my_comparisons = list(c(com[1],com[2]),
#                       c(com[3],com[4]),
#                       c(com[5], com[6]),
#                       c(com[7], com[8]),
#                       c(com[9], com[10])) 
# 
# p <- p + stat_compare_means(comparisons = my_comparisons,
#                             label = "p.signif", method = "wilcox.test", ,size = 5,paired = T)

levels(RMSE_long$variable) = gsub("_"," ",levels(RMSE_long$variable)) 
com = levels(RMSE_long$variable)
## Need to split data so single plots can be created and then combined
RMSE_tg <- RMSE_long %>%
  split(., .$variable)  

my_comparisons = list(c(com[1],com[2]),
                      c(com[3],com[4]),
                      c(com[5], com[6]),
                      c(com[7], com[8]),
                      c(com[9], com[10]))  


functionp <- function(x){ 
  p_value = x
  if (p_value <= 0.0001) {
    p_label <- "p < 0.0001"
  } else if (p_value <= 0.001) {
    p_label <- "p < 0.001"
  } else if (p_value <= 0.01) {
    p_label <- "p < 0.01"
  } else if (p_value <= 0.05) {
    p_label <- "p < 0.05"
  } else {
    p_label <- paste0("p = ", round(p_value, 3))
  } 
  p_label
}

## if needed z statistic is round(qnorm(p_values[i]/2),2)


## 4 Create plot for baseline comparisons
mplots = list()
for(i in 1:5){  
  set = my_comparisons[[i]]
  data = bind_rows(RMSE_tg[set])  
  mlabel = paste0("W = ",test_statistic[i], ", ",functionp(p_values[i]),", rB = ",round(effect_size[[i]],3))
  mplots[[i]] = ggplot(data, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(outlier.size = 1) + 
    labs(x = sub("\\ .*", "", set[1]), y = "RMSE") + 
    stat_boxplot(geom = "errorbar", width = 0.2) +
    scale_fill_discrete(name = "variable") + 
    scale_y_continuous(breaks = c(0,2,4,6,8,10,12)) +
    theme_minimal()+ 
    theme(legend.position = "none") +    
    theme(axis.text.x = element_blank()) + 
    theme(axis.title = element_text(size=16)) +
    annotate("text", x = 1.5, y = 12.5, label = mlabel)  
}  

## A separate general legend needs to created and then used as the final plot before combining
legend_plot <- ggplot() +
  geom_boxplot(aes(x = 1, y = 1, fill = "Model"), width = 0.5) +
  geom_boxplot(aes(x = 1, y = 0.7, fill = "Baseline"), width = 0.5) + 
  scale_fill_manual(values = c("Model" = "#F8766D", "Baseline" = "#00BFC4")) +
  theme_void() + theme(legend.position = c(0.5, 0.5),legend.key.size = unit(4, 'cm'), legend.title = element_blank(),legend.text = element_text(size=16),plot.background = element_rect(fill = "white")) 

mplots[[6]] = legend_plot
baseplot = patchwork::wrap_plots(mplots)     
baseplot + plot_annotation('Best model against baseline median performance', theme=theme(plot.title=element_text(hjust=0.5,size=25, face="bold"))) 
dev.print(svg,"Base_comparisons.svg",height = 10, width = 11)
dev.off()


## 5. Setup for best model comparisons
levels(RMSEcompare_long$variable) = gsub("_"," ",levels(RMSEcompare_long$variable)) 
test = kruskal.test(value ~ variable, data = RMSEcompare_long) 
hadj = unname(test$statistic)
n = sum(table(RMSEcompare_long$value,RMSEcompare_long$variable))
e2 = hadj*(n+1)/(n^2-1) ## effect size see https://peterstatistics.com/CrashCourse/3-TwoVarUnpair/NomOrd/NomOrd3c.html 

dat = RMSEcompare_long
dat$var1 = dat[,"value"]
dat$var2 = dat[,"variable"] 
dunn_stat = dat %>% dunn_test(var1 ~ var2, p.adjust.method = "fdr")
dunn_stat = dunn_stat %>% add_xy_position(x = "var2") 
dunn_stat$p.adj = sapply(dunn_stat$p.adj,functionp) 
mlabel = paste0("Kruskal-Wallis (H(4) = ",round(test$statistic,2),"), ",functionp(test$p.value),", epsilon squared = ",round(e2,3)) 


## 6. Create plot for best model comparisons
ggplot(dat, mapping = aes(x = var2, y = var1, fill = var2)) + 
  geom_boxplot() +  
  stat_boxplot(geom = "errorbar", width = 0.2) +
  stat_pvalue_manual(dunn_stat, label = "p.adj", hide.ns = F, inherit.aes = F,
                     size = 3, step.increase = 0.05) +
  labs(x = "", y = "RMSE") +  
  annotate("text", x = 2, y = 20, label = mlabel) + 
  theme_minimal() + 
  theme(legend.position = "none",axis.title = element_text(size=16),axis.text.x = element_text(size=13)) +  
  plot_annotation('Comparison of best model for each modality', theme=theme(plot.title=element_text(hjust=0.5,size=25, face="bold")))  

# dev.print(svg,"Best_comparisons.svg",height = 10, width = 11)
# dev.off()

## Different outcomes comparison 
mclin6m = final$`Sequential/Clinical`$Fsel003
mfunctional1y = final$`Sequential/FC`$Fsel021 
mmulti1y = final$`Sequential/Multimodal`$Fsel022  
mSC6m = final$`Sequential/SC`$Fsel006
mstruc6m = final$`Sequential/Structural`$Fsel001

outcomecom = list(cbind(mclin$gaussprRadialRMSE,mclin6m$gaussprRadialRMSE),cbind(mfunctional$gaussprRadialRMSE,mfunctional1y$gaussprRadialRMSE),cbind(mmulti$gaussprRadialRMSE,mmulti1y$gaussprRadialRMSE),cbind(mSC$gaussprRadialRMSE,mSC6m$gaussprRadialRMSE),cbind(mstruc$rfRMSE,mstruc6m$rfRMSE))
names(outcomecom) = c("Clinical","FC","Multimodal","SC","Morphometric") 

outcomeplots = list()
for (i in 1:5){
  data = outcomecom[[i]] 
  wilcox_result <- wilcox.test(data[,1],data[,2], paired = F) 
  p_values <- wilcox_result$p.value
  test_statistic <- wilcox_result$statistic 
  effect_size <- round(effectsize::rank_biserial(data[,1],data[,2], paired = F)[,1],3) 
  
  mlabel = paste0("U = ",test_statistic, ", ",functionp(p_values))  
  if(p_values < 0.05){mlabel = paste0(mlabel,", rB = ",effect_size)}
  data_long = melt(data)   
  data_long$Var2 = as.factor( data_long$Var2)
  levels(data_long$Var2) = c("original","opposite")
  
  outcomeplots [[i]] = ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot(outlier.size = 1) + 
    labs(x = names(outcomecom)[i], y = "RMSE") + 
    stat_boxplot(geom = "errorbar", width = 0.2) +
    scale_fill_discrete(name = "Var2") + 
    scale_y_continuous(breaks = c(0,2,4,6,8,10,12)) +
    theme_minimal()+ 
    theme(legend.position = "none") +    
    theme(axis.text.x = element_blank()) + 
    theme(axis.title = element_text(size=16)) +
    annotate("text", x = 1.5, y = 12.5, label = mlabel)  
}

legend_plot2 <- ggplot() +
  geom_boxplot(aes(x = 1, y = 1, fill = "Original\noutcome"), width = 0.5) +
  geom_boxplot(aes(x = 1, y = 0.7, fill = "Opposite\noutcome"), width = 0.5) + 
  scale_fill_manual(values = c("Original\noutcome" = "#F8766D", "Opposite\noutcome" = "#00BFC4")) +
  theme_void() + theme(legend.position = c(0.5, 0.5),legend.key.size = unit(4, 'cm'), legend.title = element_blank(),legend.text = element_text(size=16),plot.background = element_rect(fill = "white")) 


outcomeplots[[6]] = legend_plot2
baseplot = patchwork::wrap_plots(outcomeplots)     
baseplot + plot_annotation('Best model against opposite outcome model', theme=theme(plot.title=element_text(hjust=0.5,size=25, face="bold"))) 
dev.print(svg,"Outcome_comparisons.svg",height = 10, width = 11)
dev.off()

## Different functional modality comparison  
mfunctionalPASAT = final$`Sequential/FC`$Fsel003 
mmultiPASAT = final$`Sequential/Multimodal`$Fsel004

FC = cbind(mfunctional$gaussprRadialRMSE,mfunctionalPASAT$gaussprRadialRMSE)  
colnames(FC) = c("Resting-state","PASAT")
Multimodal = cbind(mmulti$gaussprRadialRMSE,mmultiPASAT$gaussprRadialRMSE) 
colnames(Multimodal) = c("Resting-state","PASAT")  

funcvar = list(FC,Multimodal) 
names(funcvar) = c("FC","Multimodal")

funcvarplots = list()
for (i in 1:2){  
  data = funcvar[[i]] 
  wilcox_result <- wilcox.test(data[,1],data[,2], paired = F) 
  p_values <- wilcox_result$p.value
  test_statistic <- wilcox_result$statistic 
  effect_size <- round(effectsize::rank_biserial(data[,1],data[,2], paired = F)[,1],3) 
  
  mlabel = paste0("U = ",test_statistic, ", ",functionp(p_values))  
  if(p_values < 0.05){mlabel = paste0(mlabel,", rB = ",effect_size)}
  data_long = melt(data)   
  data_long$Var2 = as.factor( data_long$Var2)
  
  funcvarplots [[i]] = ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot(outlier.size = 1) + 
    labs(x = names(funcvar)[i], y = "RMSE") + 
    stat_boxplot(geom = "errorbar", width = 0.2) +
    scale_fill_discrete(name = "Var2") + 
    scale_y_continuous(breaks = c(0,2,4,6,8,10,12)) +
    theme_minimal()+  
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size=7)) + 
    theme(axis.title = element_text(size=8)) +
    annotate("text", x = 1.5, y = 12.5, label = mlabel,size = 3)  
  }

baseplot = patchwork::wrap_plots(funcvarplots)     
baseplot + plot_annotation('Best model against opposite FC modality', theme=theme(plot.title=element_text(hjust=0.5,size=12.5, face="bold"))) 
dev.print(svg,"FCmodality_comparisons.svg",height = 3, width = 6)
dev.off()






