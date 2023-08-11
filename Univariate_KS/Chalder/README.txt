clean_duplicates.m: additional code that Predictors.m requires to handle duplicates from the structural connectivity results
Predictors.m: Matlab code to find predictors of a clinical variable, given neuroimaging and clinical data using general linear models (conn_glm). Outputs an Excel sheet called Results.xlsx. Includes code using change scores and controlling for baseline fatigue
Clinical.csv: table of the raw clinical data of fatigue scores from the Chalder fatigue severity questionnaire 
LIFT_data_v2.mat: final version of all neuroimaging modalities from LIFT, including structural grey matter volumes, thickness, and surface area, two functional connectivity matrices (PASAT, resting-state) and one structural connectivity matrix 
Results.xlsx: the output from Predictors.m; other Excel files identify the same results but with only PEP or CBA treatment groups, without age, gender, or imaging site as confounders, using change scores instead of controlling for baseline fatigue 
