library(dplyr)
library(tidyverse)
library(caret)   
library(DMwR2)
library(CORElearn)  
library(RGCCA)
library(kernlab)  
library(stats) 
library(reshape2)  
library(nestedcv)
library(lime) 
library(iml) 
library(yaImpute) 
library(partykit)

options(scipen = 999)
load(paste0(getwd(),"/Scipts/MLdataR_6months.RData"))
label = "Chalder2_Total" 
ccadesign = "direct" 
stability = 70 
relief = 0.25 
sigmaGrid = expand.grid(sigma = seq(1e-04, 0.015,length = 20))  
strat = "TreatmentG2" 
functional = "rs_FC1_"

preprocess_data <- function(x,y, label, ccadesign, stability, relief,functional) {
  train_cca = x
  conf = as.numeric(unclass(factor(train_cca$TBV)))
  conf = cbind(1, conf) 
  a = grep ("^Area1_17$",colnames(train_cca))  
  b = grep ("^Volume1_84$",colnames(train_cca))  
  beta_est = solve(t(conf)%*%conf)%*%t(conf)%*%as.matrix(train_cca[,a:b]) 
  train_cca[,a:b] = train_cca[,a:b] - conf%*%beta_est 
  
  ## Impute missing features using complete cases for outcome
  train_cca = DMwR2::knnImputation(train_cca)    
  
  ## Subset features according to modality
  a = grep ("^Edinburgh$",colnames(train_cca)) 
  b = grep ("^Chalder1_Total$",colnames(train_cca)) 
  x1 = train_cca[,a:b]  
  x2 = train_cca %>% dplyr::select(starts_with(functional)) 
  x3 = train_cca %>% dplyr::select(starts_with("SC1_"))  
  x4 = train_cca %>% dplyr::select(starts_with(c("Area1_","Thickness1_","Volume1_")))
  x5 = y 
  
  ## Apply RRelief for each subset using fixed nearest hit and miss neighbors (Le et al. 2019; Le, Dawkins, and McKinney 2019)   
  ## relief = 1(all positive scores), 0.25/50 - 1/4 or 1/2 of them
  myfil = function(x,y,label){
    mydata = cbind(x,y) 
    names(mydata)[length(mydata)] = label  
    k = floor((dim(mydata)[1]-1)*0.154) 
    ranked_vars = CORElearn::attrEval(label,mydata,estimator = "RReliefFequalK",
                                      outputNumericSplits=FALSE,
                                      kNearestEqual = k)
    top_vars = names(which(sort(ranked_vars, decreasing = TRUE)>0))  
    nvar = round(relief*length(top_vars))
    top_vars = top_vars[1:nvar]
    fdata = as.data.frame(x[,top_vars]) 
    colnames(fdata) = top_vars
    fdata}  
  x1 = myfil(x1,x5,label) 
  x2 = myfil(x2,x5,label) 
  x3 = myfil(x3,x5,label) 
  x4 = myfil(x4,x5,label) 
  
  ## sgCCA design: either direct (relation only with outcome) 
  ## or superblock (combine all features in one, relate to each other and then to outcome) see Garali, Imene, et al.(2018)
  if(ccadesign == "superblock"){xsup = cbind(x1,x2,x3,x4)}
  
  nsub = length(x1[,1])    
  if(ccadesign == "superblock"){ C=matrix(0,6,6);C[,5]=1;C[5, ]=1;diag(C)=0
  } else {C=matrix(0,5,5);C[1:4,5]= 1;C[5,1:4] = 1} 
  
  ## Function to calculate penalty parameter based on number of variables
  shrink = function(a){
    optimum = tau.estimate(a)
    if(is.na(optimum)){optimum = 1}
    b = length(a)
    min = 1/sqrt(b)
    if(optimum >= min){optimum} else {min}
  } 
  ## Function to extract non-zero coefficient variables
  extract = function(mycca){row.names(which(mycca!=0,arr.ind = T))} 
  
  ## Apply Stability selection: run analysis with subsamples, use features common across them
  ## Done to stabilise sgCCA results (see ng, Alex, et al.(2019)) 
  
  variables = list() 
  nsample = round(nsub*0.5)
  for(l in 1:100){
    set.seed(l)
    variableset = list()
    sample = sample(c(1:nsub), nsample, replace = F)  
    if(ccadesign == "superblock"){xlist = list(x1[sample,,drop=F],x2[sample,,drop=F],x3[sample,,drop=F],x4[sample,,drop=F],xsup[sample,,drop=F],x5[sample])
    } else {xlist = list(x1[sample,,drop=F],x2[sample,,drop=F],x3[sample,,drop=F],x4[sample,,drop=F],x5[sample])}
    if(ccadesign == "superblock"){c1 = c(sapply(xlist[-6],shrink),1)
    } else {c1 = c(sapply(xlist[-5],shrink),1)}
    cca = sgcca(xlist, C, c1,scheme = "centroid", verbose = FALSE)
    cca=cca$a[1:length(cca$a)-1] 
    if(ccadesign == "superblock"){variableset = lapply(cca[5],extract)
    } else {variableset = lapply(cca,extract)}
    if(ccadesign == "superblock"){names(variableset) = "Features" 
    } else {names(variableset) = c("Clinical","Functional","DTI","Structural")}
    variables[[l]] = variableset
  }
  if(ccadesign == "superblock"){
    variableset = lapply(variables, `[`, "Features")
    count = names(which(table(unlist(variableset)) >= stability)) 
    train_cca = as.data.frame(cbind(subset(train_cca,select = count),x5))
  } else {
    variableset1 = lapply(variables, `[`, "Clinical")
    variableset2 = lapply(variables, `[`, "Functional")
    variableset3 = lapply(variables, `[`, "DTI")
    variableset4 = lapply(variables, `[`, "Structural") 
    count1 = names(which(table(unlist(variableset1)) >= stability))
    count2 = names(which(table(unlist(variableset2)) >= stability))
    count3 = names(which(table(unlist(variableset3)) >= stability))
    count4 = names(which(table(unlist(variableset4)) >= stability)) 
    x1 = subset(x1,select = count1)
    x2 = subset(x2,select = count2)
    x3 = subset(x3, select = count3)
    x4 = subset(x4,select = count4)
    train_cca = as.data.frame(cbind(x1,x2,x3,x4))
  }
  
  ccavartemp = names(train_cca)
  # print("At ccavartemp")
  ## Add T group if not selected
  ccavar = if("TreatmentG3" %in% ccavartemp){ccavartemp} else {append(ccavartemp,"TreatmentG3",after = 0)}
  ccavar = if("TreatmentG2" %in% ccavartemp){ccavar} else {append(ccavar,"TreatmentG2",after = 0)}
  ccavar
} 

your_data = fulldata[complete.cases(fulldata$Chalder2_Total),]
y = your_data$TreatmentG2

set.seed(123)
out_folds <- caret::createFolds(y, k = 8)
in_folds <- lapply(out_folds, function(i) {
  train_y <- y[-i]
  caret::createFolds(train_y, k = 8)
})
x = your_data[,1:ncol(your_data)-1] 
y = your_data[,ncol(your_data)]   


ctrl = trainControl(method = "repeatedcv",
                    number = 8, 
                    repeats = 10,
                    search = "grid",
                    preProcOptions = c("centre","scale"))   

# res <- nestcv.train(y, x, method = "gaussprRadial",filterFUN = preprocess_data, filter_options = list(label, ccadesign, stability, relief,functional),
#                     inner_folds = in_folds,tuneGrid = sigmaGrid,
#                     outer_folds = out_folds)  

res_ctrl = nestcv.train(y, x, method = "gaussprRadial",filterFUN = preprocess_data, filter_options = list(label, ccadesign, stability, relief,functional),
                        trControl = ctrl,tuneGrid = sigmaGrid,
                        outer_folds = out_folds)  

x = x[,res_ctrl$final_vars] 

set.seed(123) 
# cors = list() 
generate_permutations <- function(data, n_permutations, max_correlation = 0.1) {
  n <- length(data)
  original_correlation <- cor(data, seq_len(n))
  permutations <- matrix(NA, nrow = n_permutations, ncol = n)
  
  for (i in 1:n_permutations) {
    perm <- sample(data)
    perm_correlation <- cor(perm, data)
    
    while (abs(perm_correlation) > max_correlation) {
      perm <- sample(data)
      perm_correlation <- cor(perm, data)
    }
    
    permutations[i,] <- perm
  }
  
  return(data.frame(permutations))
}

permuted_ys = generate_permutations(y, n_permutations = 1000, max_correlation = 0.1)
permuted_results = list() 
for (i in 1:1000) {
  # y_shuffle = sample(y)  
  # cors[[i]] = cor(y,y_shuffle) 
  y_shuffle = as.numeric(permuted_ys[i,])
  model = suppressWarnings(nestcv.train(y_shuffle, x, method = "gaussprRadial",filterFUN = NULL,
                      tuneGrid = res_ctrl$finalTune,trControl = ctrl,
                      outer_folds = out_folds,verbose = F)) 
  permuted_results[[i]] = model$summary 
}

RMSE <- sapply(permuted_results, function(element) element["RMSE"])
p_value = round(sum(RMSE < 6.74)/1000,2) 
p005 = round(RMSE[order(RMSE)][50],2)
permute_visual = as.data.frame(RMSE)
permute_visual$category <- ifelse(permute_visual$RMSE < 6.74, 'Lower error', 'Higher error') 
permute_visual$category = as.factor(permute_visual$category)

ggplot(permute_visual, aes(x = RMSE, fill = category )) +
  geom_histogram(bins = 35, color = "white") +  
  scale_fill_manual(values = c("grey","black")) +
  theme_light() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5,size=18),axis.text.x = element_text(size=12),axis.title = element_text(size=14)) + 
  labs(title = "Permuted performance distribution",fill = "Model comparison",y = "Count") + 
  geom_vline(xintercept = 6.74, linetype = "dashed", color = "blue") +  
  annotate("text", x = 7.07, y = 88, label = paste0("Model = 6.74, p = ",p_value), color = "blue",size=3) +  
  geom_vline(xintercept = p005, linetype = "dashed", color = "red") + 
  annotate("text", x = p005 - 0.3, y = 90, label = paste0("RMSE = ", p005, ", p < 0.05"), color = "red",size=3) +
  geom_vline(xintercept = mean(RMSE), linetype = "dashed", color = "black") + 
  annotate("text", x = mean(RMSE) + 0.28, y = 82, label = paste0("RMSE (mean) = ", round(mean(RMSE),2)), color = "black",size=3)

# dev.print(svg,"Permutation_test.svg",height = 5, width = 6)
# dev.off()

# calculate_median <- function(name) {
#        values <- sapply(test_acc, function(element) element[name])
#        median(values)
# }
# 
# 

# permuted_results_full = list()  
# set.seed(134) 
# for (j in 1:1000){ 
#   set.seed(j)
#   outer_folds = caret::createFolds(your_data[, strat],8, list = FALSE) 
#   permuted_data = x  
#   y_shuffle = as.numeric(permuted_ys[j,])
#   permuted_data$Chalder2_Total = y_shuffle  
#   test_acc = list()
#   for(u in 1:8){   
#     outer_idx = which(outer_folds!=u)   
#     train = permuted_data[outer_idx,]  
#     trn.pheno = train[,ncol(train)] 
#     train = train[,1:(ncol(train)-1)] 
#     test = permuted_data[-outer_idx,]  
#     tst.pheno = test[,ncol(test)] 
#     test = test[,1:(ncol(test)-1)]  
#     test = DMwR2::knnImputation(test,distData = train)
#     train = DMwR2::knnImputation(train) 
#     model = caret::train(train,trn.pheno,preProcess = c("center","scale"),method = "gaussprRadial",trControl = trainControl(method='none'),tuneGrid = res$finalTune,verbose=F) 
#     test_pred = stats::predict(model, test)  
#     test_acc[[u]] = caret::postResample(test_pred,tst.pheno)
#   } 
#   permuted_results_full[[j]] = sapply(c("RMSE", "Rsquared", "MAE"), calculate_median)
# }
# 
# RMSE_full <- sapply(permuted_results_full, function(element) element["RMSE"]) 
# p_value_full = sum(RMSE_full < 6.74)/1000  
# p005 = round(RMSE_full[order(RMSE_full)][50],2)
# permute_visual = as.data.frame(RMSE_full)
# permute_visual$category <- ifelse(permute_visual$RMSE_full < 6.723, 'Lower', 'Higher') 
# permute_visual$category = as.factor(permute_visual$category)
# 
# ggplot(permute_visual, aes(x = RMSE_full, fill = category )) +
#   geom_histogram(bins = 35, color = "white") +  
#   scale_fill_manual(values = c("grey","black")) +
#   theme_light() +
#   theme(legend.position = "bottom") + 
#   labs(title = "Permuted distribution") + 
#   geom_vline(xintercept = 6.74, linetype = "dashed", color = "blue") +  
#   annotate("text", x = 6.97 - 0.45, y = 90, label = paste0("Model = 6.74, p = ",p_value_full), color = "blue") +  
#   geom_vline(xintercept = p005, linetype = "dashed", color = "red") + 
#   annotate("text", x = p005 - 0.23, y = 90, label = paste0("RMSE = ", p005, ", p < 0.05"), color = "red") + 
#   geom_vline(xintercept = mean(RMSE_full), linetype = "dashed", color = "black") + 
#   annotate("text", x = mean(RMSE_full) + 0.2, y = 90, label = paste0("RMSE (mean) = ", round(mean(RMSE_full),2)), color = "black")

set.seed(123)
final_model = res_ctrl$final_fit
mvar = res_ctrl$final_vars

interpret_data = x  
interpret_data = DMwR2::knnImputation(x)   

predictor = Predictor$new(final_model, data = interpret_data, y = y) 
set.seed(123)
imp = FeatureImp$new(predictor, loss = "rmse",n.repetitions = 1000) 

imp$results$feature = gsub("TreatmentG2","Usual Care Group",imp$results$feature) 
imp$results$feature = gsub("TreatmentG3","CBA Group",imp$results$feature) 
imp$results$feature = gsub("Disease_duration","Disease Duration",imp$results$feature)
rois = read.csv("rois.csv",header=FALSE) 
rois = apply(rois,1,function(x){gsub("-"," ",x)}) 
toroi = function(u,z){
  text = sub("^.*?1_(.*)", "\\1",u) 
  text = unlist(strsplit(text, "_"))
  before = rois[as.integer(text[1])]
  after = rois[as.integer(text[2])]
  paste0(z," ",before," - ",after)
}
fcindex = grep("rs_FC1",imp$results$feature) 
scindex = grep("SC1",imp$results$feature)
imp$results$feature[fcindex] = sapply(imp$results$feature[fcindex],function(x){y = toroi(x,"rsFC");y})  
imp$results$feature[scindex] = sapply(imp$results$feature[scindex],function(x){y = toroi(x,"SC");y})  
imp$results$feature = gsub("Thickness1_37","Thickness LH Postcentral",imp$results$feature) 
imp$results$feature = gsub("Volume1_48","Volume LH Temporal Pole",imp$results$feature)


plot = plot(imp) 
var_importance_plot = plot + theme_light() + 
  labs(title = "Permuted feature importance of multimodal model") +  
  xlab("Feature importance (loss in RMSE)") +
  theme(plot.title = element_text(hjust = 1,size=18))   
var_importance_plot

# dev.print(svg,"Feature_importance.svg",height = 6,width = 8) 
# dev.off()

ale_depression = FeatureEffect$new(predictor, feature = "Depression")  
interact = Interaction$new(predictor)
interact_plot = plot(interact)
ale_all = FeatureEffects$new(predictor,method = 'ale')  
ale_plots = plot(ale_all) 

plotnames = ale_all$features 
plotnames = gsub("TreatmentG2","Usual Care Group",plotnames) 
plotnames = gsub("TreatmentG3","CBA Group",plotnames) 
plotnames = gsub("Disease_duration","Disease Duration",plotnames) 
fcindex = grep("rs_FC1",plotnames) 
scindex = grep("SC1",plotnames)
plotnames[fcindex] = sapply(plotnames[fcindex],function(x){y = toroi(x,"rsFC");y})  
plotnames[scindex] = sapply(plotnames[scindex],function(x){y = toroi(x,"SC");y})  
plotnames = gsub("Thickness1_37","Thickness LH Postcentral",plotnames) 
plotnames = gsub("Volume1_48","Volume LH Temporal Pole",plotnames) 

for (i in 1:length(plotnames)){
  ale_plots[[i]]$labels[[1]] = plotnames[i] 
}

list_ale = list()
for (i in 1:length(plotnames)){
  list_ale[[i]] = ale_plots[[i]]
}
names(list_ale) = plotnames 
list_ale = list_ale[imp$results$feature]

top4_plots = list_ale[c(names(list_ale)[1:4],"Usual Care Group","CBA Group")]
wrap_plots(top4_plots,ncol = 2)
# dev.print(svg,"Ale4Plots.svg",height = 6,width = 11) 
# dev.off()

other_ale = list_ale[setdiff(names(list_ale),names(top4_plots))]  
fcindex = grep("rsFC ",names(other_ale)) 
scindex = grep("SC ",names(other_ale)) 

wrap_plots(other_ale[fcindex],ncol = 2)
# dev.print(svg,"Ale_FC_Plots.svg",height = 16,width = 14) 
# dev.off() 

wrap_plots(other_ale[scindex],ncol = 2)
# dev.print(svg,"Ale_SC_Plots.svg",height = 16,width = 14) 
# dev.off() 

indices = which(!grepl("rsFC ", names(other_ale)) & !grepl("SC ", names(other_ale))) 
wrap_plots(other_ale[indices],ncol = 1)
# dev.print(svg,"Ale_other_Plots.svg",height = 4,width = 5) 
# dev.off() 

id = read.csv("ID.csv") 
id = na.omit(cbind(id,fulldata[,c("Chalder1_Total","Chalder2_Total")]))  
change = id$Chalder2_Total - id$Chalder1_Total 
id$ID[order(change)[1]]
smp_size <- floor(0.75 * nrow(x)) 
set.seed(145)
train_ind <- sample(seq_len(nrow(interpret_data)), size = smp_size)  
train_lift <- interpret_data[train_ind,]
test_lift <- interpret_data[-train_ind,] 
lift_pred <- predict(final_model, test_lift) 
change = lift_pred - your_data[-train_ind,"Chalder1_Total"] 
biggest_change <- order(change)[1]

explainer <- lime::lime(train_lift,final_model)   
set.seed(145)
explanation_single <- lime::explain(test_lift[biggest_change,], explainer,n_features = 5,n_permutations = 100000)  
explanation_single$feature = c("SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","Depression","rsFC LH Pallidum - LH Pars Orbitalis","rsFC RH Cuneus - RH Pars Triangularis","rsFC RH Caudal Anterior Cingulate Gyrus - RH Lateral Orbito Frontal Cortex") 
explanation_single$feature_desc = c("0.011 < SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","9 < Depression","rsFC LH Pallidum - LH Pars Orbitalis <= 0.019","0.26 < rsFC RH Cuneus - RH Pars Triangularis","rsFC RH Caudal Anterior Cingulate Gyrus - RH Lateral Orbito Frontal Cortex <= 0.14") 
plot = plot_features(explanation_single)  
plot$data$prediction = round(plot$data$prediction,1) 
plot$data$`Baseline` = rep(25,5) 
plot$data$`Actual follow-up` = rep(1,5) 
temp = plot[["facet"]][["params"]][["facets"]]  
temp[["Baseline"]] = quo(Baseline) 
temp[["Actual follow up"]] = quo(`Actual follow-up`)
plot[["facet"]][["params"]][["facets"]] = temp 

plot(plot)
# dev.print(svg,"Lime_single_case.svg",height = 4, width = 8) 
# dev.off()

biggest_change <- order(change)[1:5]
set.seed(145)
explanation <- lime::explain(test_lift[biggest_change,], explainer,n_features = 2,n_permutations =10000)  
explanation$feature = c("SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","Depression","SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","rsFC LH Frontal Pole - LH Insula","SC LH Medial Orbito Frontal - LH-Frontal Pole","Depression","SC LH Medial Orbito Frontal - LH-Frontal Pole","rsFC RH Cuneus - RH Pars Triangularis","rsFC LH Frontal Pole - LH Insula","Depression") 
explanation$feature_desc = c("0.11 < SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","9 < Depression","0.11 < SC RH Bank of the Superior Temporal Sulcus - RH Rostral Middle Frontal Gyrus","rsFC LH Frontal Pole - LH Insula <= -0.035","0.044 < SC LH Medial Orbito Frontal - LH-Frontal Pole","9 < Depression","0.44 < SC LH Medial Orbito Frontal - LH-Frontal Pole","rsFC RH Cuneus - RH Pars Triangularis  <= -0.067","0.184 < rsFC LH Frontal Pole - LH Insula","9 < Depression") 
plot_explanations(explanation)  
# dev.print(svg,"Lime_5_cases.svg",height = 4, width = 8) 
# dev.off()


