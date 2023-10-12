.libPaths(c("/uoa/scratch/localscratch/users/rdo071/LIFT/Kris_Cluster/R_Library",.libPaths())) 

# 1. Libraries

library(tidyverse)
library(doParallel)
library(caret)   
library(glmnet) 
library(randomForest)  
library(gbm)
library(xgboost)
library(DMwR2)
library(CORElearn)  
library(RGCCA)
library(kernlab)  
library(modelgrid) 
library(stats) 
library(mice)
library(reshape2) 
library(doSNOW)

# 2. Algorithms, parameter settings. feature selecion, and CV

## Set parameters for script
label = "Chalder2_Total" ## outcome
strat = "TreatmentG2" ## used to stratify to folds
algor = c("glmnet","rf","svmLinear","svmRadial","gaussprRadial","xgbDART","gbm","Base") ## algorithm names
## Hyper-parameter grids
svmGrid = expand.grid(C = c(1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1,1,10,100,1000,2000))  
sigmaGrid = expand.grid(sigma = seq(1e-04, 0.015,length = 20))  
## Settings of feature selection:
## Relief - percentage of positive ranked v 
## Stability - percentage of subsamples that have variables (for CCA) 
## Canonical correlation analysis (CCA) structure - direct or super-block 
f_sel = expand.grid(relief = c(0.25,0.5,1),stability = c(50,70,90),fdata =c("PASAT","REST"),ddata = c("6m","1y"))
## Cross-validation parameters
nrepeats = 10  
ncv_folds = c(10,10)  
nimputed = 10 ## for multiple imputation  
modelsave = list()

# 3. ML code

# ptm = proc.time() 

finalf = list() 
finalf = for(u in 1:nrow(f_sel)) {  
  relief = f_sel[u,]$relief 
  stability = f_sel[u,]$stability  
  fdata <- f_sel[u,]$fdata   
  ddata = f_sel[u,]$ddata
  if (fdata == "PASAT" & ddata == "6m"){load(paste0(getwd(),"/MLdataP_6months.RData"));functional = "P_FC1_"} else if (fdata == "REST" & ddata == "6m"){load(paste0(getwd(),"/MLdataR_6months.RData"));functional = "rs_FC1_"} else if (fdata == "PASAT" & ddata == "1y"){load(paste0(getwd(),"/MLdataP_1year.RData"));functional = "P_FC1_"} else if(fdata == "REST" & ddata == "1y") {load(paste0(getwd(),"/MLdataR_1year.RData"));functional = "rs_FC1_"}
  
  ## To regress Total Brain Volume (TBV) from structural measures
  conf = as.numeric(unclass(factor(fulldata$TBV))) 
  conf = cbind(1, conf)  
  
  ## Create data frame with performance metrics for each algorithm 
  methods = algor
  suffix = c("RMSE","R2","MAE")  
  names = paste(rep(methods, each = length(suffix)), suffix, sep = "") 
  all_results = as.data.frame(matrix(nrow=1, ncol = length(names))) 
  names(all_results) = names
  for(z in 1:nrepeats){
    set.seed(z) 
    results = as.data.frame(matrix(nrow=1, ncol = length(names)))  
    names(results) = names
    outer_folds = caret::createFolds(fulldata[, strat], ncv_folds[1], list = FALSE)   
    
    for (i in 1:ncv_folds[1]){  
      inner_folds = caret::createFolds(fulldata[, strat][outer_folds!=i], ncv_folds[2], list = TRUE)
      ## Parallel for feature selection and hyperparameter tunning
      cl=makeCluster(10) 
      registerDoParallel(cl)  
      atts <- foreach(j=1:ncv_folds[2],.packages = c('foreach', 'caret',"DMwR2","dplyr",'CORElearn','RGCCA',"kernlab","gbm","modelgrid","stats","mice",'tidyverse')) %dopar% {  
        inner_idx = which(outer_folds!=i)[-inner_folds[[j]]]   
        test_idx = which(outer_folds!=i)[inner_folds[[j]]]  
        train_cca = fulldata[inner_idx,]  
        
        ## Regress TBV for each inner fold
        C_train = conf[inner_idx, ]
        C_test = conf[test_idx, ]
        a = grep ("^Area1_17$",colnames(train_cca))  
        b = grep ("^Volume1_84$",colnames(train_cca))  
        beta_est = solve(t(C_train)%*%C_train)%*%t(C_train)%*%as.matrix(train_cca[,a:b]) 
        train_cca[,a:b] = train_cca[,a:b] - C_train%*%beta_est 
        
        train_cca = train_cca[complete.cases(train_cca[,label]),] 
        
        ## Subset features to Structural and outcome
        x1 = train_cca %>% dplyr::select(starts_with(c("Area1_","Thickness1_","Volume1_")))
        x5 = train_cca[,label] 
        
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
        
        ## sgCCA design: direct (relation only with outcome) 
        
        nsub = length(x1[,1])   
        C=matrix(0,2,2);C[1:1,2]= 1;C[2,1:1] = 1 
        
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
          xlist = list(x1[sample,,drop=F],x5[sample])
          c1 = c(shrink(xlist[[1]]),1)
          cca = sgcca(xlist, C, c1,scheme = "centroid", verbose = FALSE)
          cca=cca$a[1:length(cca$a)-1] 
          variableset = lapply(cca,extract)
          names(variableset) = "Structural"
          variables[[l]] = variableset
        }
        
        variableset1 = lapply(variables, `[`, "Structural")
        count1 = names(which(table(unlist(variableset1)) >= stability))
        x1 = subset(x1,select = count1)
        train_cca = as.data.frame(cbind(x1,x5))
        names(train_cca)[length(train_cca)] = label 
        ccavartemp = names(train_cca)  
        ## Add T group if not selected
        ccavar = if("TreatmentG3" %in% ccavartemp){ccavartemp} else {append(ccavartemp,"TreatmentG3",after = 0)}
        ccavar = if("TreatmentG2" %in% ccavartemp){ccavar} else {append(ccavar,"TreatmentG2",after = 0)}  
        
        ## Regress TBV from structural in train, use model in test
        mltrain = fulldata[inner_idx,] 
        mltest = fulldata[test_idx,] 
        a = grep ("^Area1_17$",colnames(mltrain))  
        b = grep ("^Volume1_84$",colnames(mltrain))   
        mltrain[,a:b] = mltrain[,a:b] - C_train%*%beta_est 
        mltest[,a:b] = mltest[,a:b] - C_test%*%beta_est
        ## Create train and test split from complete outcome cases, extracte outcome
        mltrain = mltrain[complete.cases(mltrain[,label]),ccavar]  
        mltest = mltest[complete.cases(mltest[,label]),ccavar]  
        trn.pheno = mltrain[,label] 
        tst.pheno = mltest[,label]     
        mltrain = subset(mltrain, select = -grep(label, colnames(mltrain)))
        mltest = subset(mltest, select = -grep(label, colnames(mltest)))  
        trn.data = as.matrix(mltrain)
        tst.data = as.matrix(mltest)  
        
        ## Run hyperparameter tunning using same rCV, apply to test set and select optimal hyperparameters and features
        set.seed(123)
        seeds = list(vector(mode = "list", length = 51)) 
        for(s in 1:50){seeds[[s]] = sample.int(4000, 3600)} 
        seeds[[51]] = sample.int(4000, 1)  
        myindex = createMultiFolds(trn.data[,strat], 10, 5)
        mg = 
          model_grid() %>%
          share_settings(
            y = trn.pheno,
            x = trn.data,
            metric = "RMSE",
            trControl = trainControl(preProc = list(c("center","scale")),
                                     method = "repeatedcv", 
                                     number = 10, repeats = 10, 
                                     search = "grid",index = myindex,seeds = seeds,allowParallel = F,predictionBounds = c(0,33))) 
        
        model_structure = mg %>%
          add_model(model_name = "glmnet",
                    method = "glmnet",tuneLength = 10) %>%
          add_model(model_name = "rf", 
                    method = "rf",tuneGrid = expand.grid(mtry = round(seq(1,floor(ncol(trn.data)/3),length = 10)))) %>%
          add_model(model_name = "svmLinear",
                    method = "svmLinear",tuneGrid = svmGrid) %>%  
          add_model(model_name = "svmRadial",
                    method = "svmRadial",tuneLength = 10) %>%  
          add_model(model_name = "gaussprRadial",
                    method = "gaussprRadial",tuneGrid = sigmaGrid) %>% 
          add_model(model_name = "xgbDART",
                    method = "xgbDART",tuneLength = 1,objective="reg:squarederror",nthread = 1) %>%
          add_model(model_name = "gbm",
                    method = "gbm", tuneLength = 10,verbose = F)
        model_list = train(model_structure,resample_seed = 123)  
        
        p = as.data.frame(model_list$model_fits %>% predict(., newdata = tst.data))  
        Accuracy = list()
        for (n in 1:ncol(p)){
          Accuracy[[n]] = caret::postResample(p[[n]],tst.pheno)[[1]]
        } 
        names(Accuracy) = colnames(p)   
        accu = Accuracy$glmnet
        glmnet = rbind(data.frame(accu,model_list$model_fits$glmnet$bestTune)) 
        accu = Accuracy$rf
        rf = rbind(data.frame(accu,model_list$model_fits$rf$bestTune))  
        accu = Accuracy$svmLinear
        svmLinear = rbind(data.frame(accu,model_list$model_fits$svmLinear$bestTune)) 
        accu = Accuracy$svmRadial
        svmRadial = rbind(data.frame(accu,model_list$model_fits$svmRadial$bestTune)) 
        accu = Accuracy$gaussprRadial
        gaussprRadial = rbind(data.frame(accu,model_list$model_fits$gaussprRadial$bestTune))   
        accu = Accuracy$xgbDART
        xgbDART = rbind(data.frame(accu,model_list$model_fits$xgbDART$bestTune))
        accu = Accuracy$gbm
        gbm = rbind(data.frame(accu,model_list$model_fits$gbm$bestTune))
        mylist = list(ccavar,glmnet,rf,svmLinear,svmRadial,gaussprRadial,xgbDART, gbm) 
        methods = c("glmnet","rf","svmLinear","svmRadial","gaussprRadial","xgbDART","gbm")  
        listnames = c("Variables",methods)
        names(mylist) = listnames 
        mylist
      }
      stopCluster(cl)
      registerDoSEQ()
      ## Function to extract best hyperparameters and features for each model
      model_sel = function(a){
        s_variables = lapply(atts, `[`, "Variables") 
        params = lapply(atts, `[`, a) 
        tempvar = rbind.data.frame(do.call(rbind, s_variables)) 
        tempparam = do.call(rbind,do.call(rbind, params)) 
        vpar = cbind(tempparam, tempvar)
        best_vpar = vpar[which.min(vpar$accu), ] 
        best_var = unlist(best_vpar[,"Variables"],use.names=FALSE) 
        selected_v = if("TreatmentG3" %in% best_var){best_var} else {append(best_var,"TreatmentG3",after = 0)}
        selected_v = if("TreatmentG2" %in% best_var){selected_v} else {append(selected_v,"TreatmentG2",after = 0)} 
        best_vpar$Variables = list(c(selected_v)) 
        best_vpar
      }
      
      models = sapply(methods[1:7],model_sel)
      
      outer_idx = which(outer_folds!=i)   
      train_mod = fulldata[outer_idx,]
      test_mod = fulldata[-outer_idx,]  
      ## Regress TBV for outer fold
      C_train = conf[outer_idx, ]
      C_test = conf[-outer_idx, ]
      a = grep ("^Area1_17$",colnames(train_mod))  
      b = grep ("^Volume1_84$",colnames(train_mod))  
      beta_est = solve(t(C_train)%*%C_train)%*%t(C_train)%*%as.matrix(train_mod[,a:b]) 
      train_mod[,a:b] = train_mod[,a:b] - C_train%*%beta_est 
      test_mod[,a:b] = test_mod[,a:b] - C_test%*%beta_est   
      
      
      set.seed(321) 
      ## Create training set based on best features, apply multiple imputation
      ## and then train model using best hyperparameters; test model using imputed data sets
      apply_model = function(Variables,hyparam,smodel){
        train = train_mod[,unlist(Variables)] 
        test = test_mod[unlist(Variables)] 
        train$TreatmentG2 = as.factor(train$TreatmentG2)
        train$TreatmentG3 = as.factor(train$TreatmentG3)
        test$TreatmentG2 = as.factor(test$TreatmentG2)
        test$TreatmentG3 = as.factor(test$TreatmentG3) 
        m = nimputed  
        emptyim = mice(train,maxit = 0)   
        post = emptyim$post
        post["Chalder2_Total"] ="imp[[j]][, i] = squeeze(imp[[j]][, i], c(0, 33))" 
        imp = mice(data = rbind(train,test), ignore = c(rep(FALSE, nrow(train)), rep(TRUE,nrow(test))),maxit = 20,m = m,seed = 123,print = F,post= post,method = "rf") 
        fullimp = list()
        train_imp = list() 
        test_imp = list() 
        for (k in seq(1:imp$m))
        {
          fullimp[[k]] = complete(imp,k) 
          indx = sapply(fullimp[[k]], is.factor) 
          fullimp[[k]][indx] = lapply(fullimp[[k]][indx], function(x) as.numeric(as.character(x))) 
          train_imp[[k]] = fullimp[[k]][1:nrow(train),] 
          test_imp[[k]] = fullimp[[k]][-(1:nrow(train)),] 
        } 
        myimputations = list(train_imp,test_imp) 
        names(myimputations) = c("train","test")   
        
        accuraries = data.frame(matrix(ncol=3,nrow=1))
        for (t in 1:m){
          trn.pheno = myimputations$train[[t]][,label]   
          tst.pheno = myimputations$test[[t]][,label] 
          trn.data = subset(myimputations$train[[t]], select = -grep(label, colnames(myimputations$train[[t]])))
          trn.data = as.matrix(trn.data) 
          tst.data = subset(myimputations$test[[t]], select = -grep(label, colnames(myimputations$test[[t]])))
          trn.data = as.matrix(trn.data)    
          
          model = caret::train(trn.data,trn.pheno,preProcess = c("center","scale"),method = smodel,trControl = trainControl(method='none'),tuneGrid = hyparam,verbose=F) 
          test_pred = stats::predict(model, tst.data) 
          test_acc = caret::postResample(test_pred,tst.pheno)
          accuraries[t,] = test_acc
        }
        accuraries
      }
      
      cl=makeCluster(length(models)) 
      registerDoParallel(cl)  
      modelper <- foreach(q=1:length(models),.combine = cbind,.packages = c('foreach', 'caret',"DMwR2","dplyr","kernlab","gbm","stats","mice",'tidyverse'),.export = c("apply_model")) %dopar% {   
        var = models[[q]]$Variables 
        hyp = length(models[[q]]) - 1
        mmodel = apply_model(var,models[[q]][2:hyp],methods[q])
      }
      stopCluster(cl)
      registerDoSEQ() 
      
      mltrain = train_mod[,c("Chalder1_Total",label)] 
      mltest = test_mod[,c("Chalder1_Total",label)]
      
      m = nimputed
      emptyim = mice(mltrain,maxit = 0)  
      post = emptyim$post
      post["Chalder2_Total"] ="imp[[j]][, i] = squeeze(imp[[j]][, i], c(0, 33))"
      imp = mice(rbind(mltrain,mltest), ignore = c(rep(FALSE, nrow(mltrain)), rep(TRUE,nrow(mltest))),maxit = 20,m = m,seed = 123,print = F,post= post,method = "rf") 
      fullimp = list()
      train_imp = list() 
      test_imp = list()
      for (k in seq(1:imp$m))
      {
        fullimp[[k]] = complete(imp,k) 
        train_imp[[k]] = fullimp[[k]][1:nrow(mltrain),] 
        test_imp[[k]] = fullimp[[k]][-(1:nrow(mltrain)),] 
      } 
      myimputations = list(train_imp,test_imp) 
      names(myimputations) = c("train","test")   
      Base = data.frame(matrix(ncol=3,nrow=1))
      for (t in 1:m){
        trn.pheno = myimputations$train[[t]][,label]    
        tst.pheno = myimputations$test[[t]][,label]  
        test_pred_base = rep(median(trn.pheno), each = length(tst.pheno)) 
        test_acc_base = caret::postResample(test_pred_base,tst.pheno)  
        Base[t,] = test_acc_base
      }
      
      foldresults = cbind(modelper,Base) 
      colnames(foldresults) = colnames(results)   
      tempsave = models 
      tempsave[[8]] = foldresults
      names(tempsave) = c("glmnet","rf","svmLinear","svmRadial","gaussprRadial","xgbDART","gbm","results")
      save(foldresults,file = paste0(getwd(),"/Saved/",paste0("Fsel",u,"_Repeat",z,"_Fold",i,".RData")))
      modelsave = append(modelsave,list(tempsave)) 
      names(modelsave)[length(modelsave)] = paste0("Fsel",u,"_Repeat",z,"_Fold",i)
      results = rbind(results,foldresults)
    } 
    results = results [-1,]
    all_results = rbind(all_results,results)  
  } 
  all_results = all_results [-1,] 
  paste(paste0("Fsel_",u),"is done out of",nrow(f_sel))
  finalf[[u]] <- all_results
} 

# elapsed.time = (proc.time() - ptm)[3] 
# duration = elapsed.time/3600 
# duration

# 4. Plots

for (v in 1:nrow(f_sel)){names(finalf)[v] = paste(f_sel[v,]$relief,f_sel[v,]$stability,f_sel[v,]$ccadesign,f_sel[v,]$fdata,f_sel[v,]$ddata)}
methods = algor
suffix = c("Median","IQR")
names = paste(rep(methods, each = length(suffix)), suffix, sep = "")
summary = as.data.frame(matrix(nrow=nrow(f_sel), ncol = length(names)))
names(summary) = names
for (v in 1:nrow(f_sel)){rownames(summary)[v] = paste(f_sel[v,]$relief,f_sel[v,]$stability,f_sel[v,]$ccadesign,f_sel[v,]$fdata,f_sel[v,]$ddata)}

mindex = seq(1,ncol(finalf[[1]]),3)
cindex = seq(1,ncol(summary),2)

for (v in 1:length(mindex)){
  n = cindex[v]
  summary[,n] = apply(sapply(finalf,"[",,mindex[v]),2,median)
  summary[,n+1] = apply(sapply(finalf,"[",,mindex[v]),2,IQR)
}


finalf2 = finalf  
for (v in 1:nrow(f_sel)){
  final = finalf2[[v]] 
  combined = final %>% select(-contains(c("MAE","R2")))  
  finalf2[[v]] = combined
}


combined = melt(finalf2,)
names(combined) = c("Algorithm","RMSE","Fselect")
levels(combined$Algorithm) = algor

median_cl_boot = function(x, conf = 0.95) {
  lconf = (1 - conf)/2
  uconf = 1 - lconf
  require(boot)
  bmedian = function(x, ind) median(x[ind])
  bt = boot(x, bmedian, 10000)
  bb = boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t,
                                                                          uconf))
}

plot = ggplot(combined, aes(x=Algorithm, y=RMSE, group=Algorithm)) +
  geom_boxplot(aes(color=Algorithm)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  theme_bw()+ facet_wrap("Fselect")+
  stat_summary(fun.data = median_cl_boot, geom = "errorbar")+
  stat_summary(fun = median, geom = "point") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title = element_text(color="black", size=24, face="bold"),
        strip.text.x = element_text(color="black", size=15, face="plain"),
        axis.text = element_text(color="black", size=12, face="plain"),
        strip.background = element_rect(colour="black", fill="white",
                                        size = 1.5,linetype="solid"),
        plot.subtitle=element_text(size=5, hjust=0.5, face="plain",
                                   color = "black"))  
finalsave = list(modelsave,finalf,summary,plot) 
names(finalsave) = c("Models","Performance","Summary","plot")
save(finalsave,file = paste0(getwd(),"/Sequential_Structural.RData"))

