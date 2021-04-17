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

%%% Extract treatment groups 
% tgroups = double(~group(:,2)).';  %%% create joined treatment groups
data = data([data.group] ~= 2); %%% exclude group  
group = extractfield(data,'group');
group = dummyvar(group);  
group1 = group(:,1).'; 
group2 = group(:,2).'; 
group3 = group (:,3).';

%%% Load flag (for indexing missing data)
flag=extractfield(data,'flag');
flag=reshape(flag,10,[]);

%%% Extract covariates data
gender=extractfield(data,'gender'); 
TIV = str2double(string(extractfield(data,'TIV')));  
age=extractfield(data,'age');
site=extractfield(data,'site');
chldr1=extractfield(data,'chldr1');
chldr2=extractfield(data,'chldr2'); 
chldrdif= chldr2 - chldr1; 
pain1 = extractfield(data,'Pain1');
pain2 = extractfield(data,'Pain2'); 
%%% Create mask to include subjects with non-missing data
% Mask for Pain, first line for follow-up pain, second for Conn
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
%%% Group 2(control), change the other groups for different contrasts
% followUp1 = pain2(mask)'.*group1(mask)';
% followUp2 = pain2(mask)'.*group3(mask)'; 
followUp1 = chldr2(mask)'.*group1(mask)'; 
followUp2 = chldr2(mask)'.*group3(mask)';
covName=['control' 'treatment' 'gender' 'age' 'site' 'TIV' 'baseline' 'control_followUp' 'treatment_followUp']; 
% X=[group1(mask)',group3(mask)', gender(mask)', age(mask)',
% site(mask)',TIV(mask)',chldr1(mask)',followUp1,followUp2]; 
X=[ gender(mask)', age(mask)', site(mask)',TIV(mask)',followUp1,followUp2];
% X=[group1(mask)',group3(mask)', gender(mask)', age(mask)', site(mask)',TIV(mask)',pain1(mask)',followUp1,followUp2]; 
% C=[0 0 0 0 0 0 0 -1 1]; %%% to evaluate whether the follow-up fatigue
% effect is different between the two groups    
C=[0 0 0 0 -1 1];


%%% to create GLM using the covariate and Volume data
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