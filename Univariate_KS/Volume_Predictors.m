clc
clear

%%% Load neuroimaging and covariate data (including Chalder) and Pain
%%% Severity data seperately 
load('C:\Users\Neuroimflammation-2\Kris\OneDrive - University of Glasgow\Salim\LIFT_data_v2.mat'); 
Pain = readtable('C:\Users\Neuroimflammation-2\Kris\OneDrive - University of Glasgow\Salim\LIFT\Imputations\2020-05-28_LIFT_clinical\Pain.csv'); 
for n=1:88 
    data(n).Pain1 = Pain.Pain1(n); 
    data(n).Pain2 = Pain.Pain2(n);  
end

%%% Extract group 
data = data([data.group] ~= 2); %%% exclude group
% data = data([data.group] == 3); %%% include group1 only 
%%% Load flag (for indexing missing data)
flag=extractfield(data,'flag');
flag=reshape(flag,10,[]);

%%% Extract covariates data
gender=extractfield(data,'gender');  
TIV = str2double(string(extractfield(data,'TIV')));  
TBV = extractfield(data,'TBV');
age=extractfield(data,'age');
site=extractfield(data,'site');
chldr1=extractfield(data,'chldr1');
chldr2=extractfield(data,'chldr2'); 
chldrdif= chldr2 - chldr1; 
pain1 = extractfield(data,'Pain1');
pain2 = extractfield(data,'Pain2'); 
paindif = pain2 - pain1;
%%% Create mask to include subjects with non-missing data 
% Mask for Chalder/Pain, first line for follow-up chalder, second for
% Volume
% mask = ~isnan(pain2); 
mask = ~isnan(chldr2); 
mask = flag(9,:) + mask==2;  

%%% Extract connectivity data
SC1=[];
subj=1;
name=[]; 
for i=1:length(data)
    if mask(i)==1
        name(subj)=(data(i).name);
        SC1(:,:,subj)=data(i).Volume1;
        subj=subj+1;
    end
end

%%%%%%% GLM construction %%%%%%

%%% Design and contrast construction (all scenarios)
covName=['subj' 'gender' 'age' 'site' 'TIV' 'outcome']; %%% unadjusted/difference
% covName=['subj' 'gender' 'age' 'site' 'baseline' 'outcome']; %%% adjusted
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',chldr2(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',chldrdif(mask)'];
X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',chldr1(mask)',chldr2(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',pain2(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',paindif(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',TBV(mask)',pain1(mask)',pain2(mask)'];
% C=[0 0 0 0 0 1]; 
C=[0 0 0 0 0 0 1]; 

%%% to create GLM using the covariate and Conn data
LME.xX.X=X(setdiff(1:(subj-1),0),:);
LME.xX.name=covName;
LME.y=squeeze(SC1(:,:,setdiff(1:(subj-1),0)))';
LME.c=C;
[h,F,p,dof, statsname, B]=conn_glm(LME.xX.X, LME.y, LME.c, [1]);
p=2*min(p,1-p); %%% create two-tailed p value
LME.p_FDR=conn_fdr(squeeze(p));    
LME.Effect_size = h; 
LME.Test_statistic = F; 
LME.Raw_Pvalue = p;
        
%%% If loop to index connections below threshold 
for i=1:length(LME.p_FDR)
    if LME.p_FDR(i) <0.05
        signif{i} = 'YES';
    else
        signif{i} = 'NO'; 
    end
end   
%%% Index for significant connections
Index = find(contains(signif,'YES')); 

%%% Create table to run stepwise regression and plot 
% seed=11;
% Gender = gender(mask)'; 
% Age = age(mask)';  
% Site = site(mask)'; 
% Chalder1 = chldr1(mask)'; 
% Chalder2 = chldr2(mask)';
% Chalderdif = chldrdif(mask)'; 
% cTIV = TIV(mask)';
% Vol = LME.y(:,seed);  
% % tbl = table(Gender,Age,Site,cTIV,Vol,Chalderdif,'VariableNames',{'Gender','Age','Site','TIV','Vol','Chalderdif'}); 
% % lowerMdl = 'Chalderdif ~ 1 + Vol';
% % md = stepwiselm(tbl,lowerMdl,'ResponseVar','Chalderdif','PredictorVars',{'Gender','Age','Site','TIV','Vol'},...
% % 'CategoricalVar',{'Gender','Site'},'Verbose',2,'PRemove',1)   
% tbl = table(Gender,Age,Site,cTIV,Vol,Chalder1,Chalder2,'VariableNames',{'Gender','Age','Site','TIV','Vol','Chalder1','Chalder2'}); 
% lowerMdl = 'Chalder2 ~ 1 + Vol';
% md = stepwiselm(tbl,lowerMdl,'ResponseVar','Chalder2','PredictorVars',{'Gender','Age','Site','TIV','Vol','Chalder1'},...
% 'CategoricalVar',{'Gender','Site'},'Verbose',2,'PRemove',1)   
