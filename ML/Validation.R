library(tidyverse)
library(caret)   
library(DMwR2)
library(CORElearn)  
library(RGCCA)
library(kernlab)  
library(stats) 
library(reshape2)  
library(nestedcv)
library(yaImpute) 
library(partykit)

load("All_scripts.RData")  
f_sel = expand.grid(relief = c(0.25,0.5,1),stability = c(50,70,90),ccadesign = c("direct","superblock"),fdata =c("PASAT","REST"),ddata = c("6m","1y"))
f_sel = filter(f_sel, ccadesign != "superblock" | stability != 50 & stability != 70 | relief!=0.5)
f_sel = filter(f_sel,stability != 90 & ccadesign != "superblock") 
row.names(f_sel) = paste0("Fsel",sprintf('%0.3d',as.numeric(row.names(f_sel))))
f_sel = filter(f_sel,fdata == "PASAT" & ddata == "6m")  
PASAT = final$`Sequential/Multimodal`[row.names(f_sel)] 
PASAT_median = data.frame(sapply(PASAT, function(x){apply(x,2,median)}))  

PASAT_RMSE = PASAT_median[grepl("RMSE",row.names(PASAT_median)),] 
PASAT_RMSE = PASAT_RMSE[1:7,]
RMSE_value <- min(PASAT_RMSE,na.rm =T) 
min_indices <- which(PASAT_RMSE == RMSE_value, arr.ind = TRUE)
RMSE_algorithm <- gsub("RMSE","",row.names(PASAT_RMSE)[min_indices[1, 1]])
RMSE_f_sel <- colnames(PASAT_RMSE)[min_indices[1, 2]]

PASAT_R2 = PASAT_median[grepl("R2",row.names(PASAT_median)),] 
PASAT_R2 = PASAT_R2[1:7,]
R2_value <- max(PASAT_R2,na.rm =T) 
max_indices <- which(PASAT_R2 == R2_value, arr.ind = TRUE)
R2_algorithm <- gsub("R2","",row.names(PASAT_R2)[max_indices[1, 1]])
R2_f_sel <- colnames(PASAT_R2)[max_indices[1, 2]] 

PASAT_MAE = PASAT_median[grepl("MAE",row.names(PASAT_median)),] 
PASAT_MAE = PASAT_MAE[1:7,]
MAE_value <- min(PASAT_MAE,na.rm =T) 
min_indices <- which(PASAT_MAE == MAE_value, arr.ind = TRUE)
MAE_algorithm <- gsub("MAE","",row.names(PASAT_MAE)[min_indices[1, 1]])
MAE_f_sel <- colnames(PASAT_MAE)[min_indices[1, 2]]


load(paste0(getwd(),"/Scipts/MLdataP_6months.RData"))
label = "Chalder2_Total" 
ccadesign = "direct" 
stability = f_sel[RMSE_f_sel,2] 
relief = f_sel[RMSE_f_sel,1] 
sigmaGrid = expand.grid(sigma = seq(1e-04, 0.015,length = 20))  
strat = "TreatmentG2" 
functional = "P_FC1_" 
algorithm = RMSE_algorithm

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
y = your_data$Chalder2_Total
x = your_data[,1:(ncol(your_data)-1)]
features = preprocess_data(x,y,label,ccadesign,stability,relief,functional) 

y = your_data$TreatmentG2
out_folds <- caret::createFolds(y, k = 8)
y = your_data$Chalder2_Total 

ctrl = trainControl(method = "repeatedcv",
                    number = 8, 
                    repeats = 10,
                    search = "grid",
                    preProcOptions = c("centre","scale"))   

res_ctrl = nestcv.train(y, x, method = RMSE_algorithm,filterFUN = preprocess_data, filter_options = list(label, ccadesign, stability, relief,functional),
                        trControl = ctrl,tuneGrid = sigmaGrid,
                        outer_folds = out_folds)  

final_model_data = x[,features]   
final_model_data[,3:ncol(final_model_data)] = scale(final_model_data[,3:ncol(final_model_data)])
final_model_data = cbind(final_model_data,y) 
colnames(final_model_data)[ncol(final_model_data)] = "Chalder2_Total"  


set.seed(123)
final_model <- train(Chalder2_Total~ .,data=final_model_data,
  method = RMSE_algorithm, 
  tuneGrid = res_ctrl$finalTune,
  trControl = trainControl(method = 'none'),   
  metric = "RMSE")  

load(paste0(getwd(),"/Scipts/MLdataP_6months_Pfizer.RData"))   
pfizer = fulldata[,features]  
pfizer$Chalder2_Total = fulldata$Chalder2_Total  
pfizer = na.omit(pfizer) 
y = pfizer$Chalder2_Total 
pfizer = pfizer[,features]   
pfizer[,3:ncol(pfizer)] = scale(pfizer[,3:ncol(pfizer)]) 
pfizer = cbind(pfizer,y)
colnames(pfizer)[ncol(pfizer)] = "Chalder2_Total"  

fitpfizer <- predict(final_model, newdata = pfizer) 
performance_pfizer = caret::postResample(fitpfizer,pfizer$Chalder2_Total) 
