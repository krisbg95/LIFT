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
% data = data([data.group] ~= 2); %%% exclude control group
data = data([data.group] == 3); %%% include one group only 
%%% Load flag (for indexing missing data)
flag=extractfield(data,'flag');
flag=reshape(flag,10,[]);

%%% Extract covariates data
gender=extractfield(data,'gender'); 
age=extractfield(data,'age');
site=extractfield(data,'site');
chldr1=extractfield(data,'chldr1');
chldr2=extractfield(data,'chldr2'); 
chldrdif= chldr2 - chldr1; 
pain1 = extractfield(data,'Pain1');
pain2 = extractfield(data,'Pain2'); 
paindif = pain2 - pain1;
%%% Create mask to include subjects with non-missing data
%%% 1(Chalder1), 2(Chalder2), 3(P_FC1)...
% Mask for Pain/Chalder, first line for follow-up, second for Conn
% mask = ~isnan(pain2); 
mask = ~isnan(chldr2); 
mask = flag(3,:) + mask==2;  

%%% Extract connectivity data
SC1=[];
subj=1;
name=[]; 
for i=1:length(data)
    if mask(i)==1
        name(subj)=(data(i).name);
        SC1(:,:,subj)=data(i).P_FC1;
        subj=subj+1;
    end
end

%%%%%%% GLM construction %%%%%%

%%% Design and contrast construction (all scenarios)
% covName=['subj' 'gender' 'age' 'site' 'outcome']; %%% unadjusted/difference
covName=['subj' 'gender' 'age' 'site' 'baseline' 'outcome']; %%% adjusted
X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',chldr2(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',chldrdif(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',chldr1(mask)',chldr2(mask)'];

% X=[ones(subj-1,1), gender(mask)',chldr1(mask)',chldr2(mask)'];
% X=[ones(subj-1,1), age(mask)',chldr1(mask)',chldr2(mask)'];  
% X=[ones(subj-1,1), site(mask)',chldr1(mask)',chldr2(mask)'];

% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',pain2(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)', site(mask)',paindif(mask)'];
% X=[ones(subj-1,1), gender(mask)', age(mask)',site(mask)',pain1(mask)',pain2(mask)'];  

% X=[ones(subj-1,1), gender(mask)',pain1(mask)',pain2(mask)'];
% X=[ones(subj-1,1), age(mask)',pain1(mask)',pain2(mask)'];  
% X=[ones(subj-1,1), site(mask)',pain1(mask)',pain2(mask)'];


% C=[0 0 0 1]; 
C=[0 0 0 0 1]; 
% C=[0 0 0 0 0 1]; 

%%% For loop through each ROI using all subjects 
%%% to create GLM using the covariate and Conn data
for ROI=1:84
        LME(ROI).xX.X=X(setdiff(1:(subj-1),0),:);
        LME(ROI).xX.name=covName;
        LME(ROI).y=squeeze(SC1(ROI,:,setdiff(1:(subj-1),0)))';
        LME(ROI).c=C;
        [h,F,p,dof, statsname, B]=conn_glm(LME(ROI).xX.X, LME(ROI).y, LME(ROI).c, [1]);
        p=2*min(p,1-p); %%% create two-tailed p value
        LME(ROI).p_FDR=conn_fdr(squeeze(p));    
        LME(ROI).Effect_size = h; 
        LME(ROI).Test_statistic = F; 
        LME(ROI).Raw_Pvalue = p;
        %%% If loop to index connections below threshold 
        if sum(LME(ROI).p_FDR <0.05)>0 
            signif{ROI} = 'YES'; 
        else 
            signif{ROI} = 'NO';
        end   
end
%%% Index for significant connections
Index = find(contains(signif,'YES')); 
