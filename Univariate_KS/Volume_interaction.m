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

%%% Extract groups 
group = extractfield(data,'group');
group = dummyvar(group);  
% tgroups = double(~group(:,2)).'; 
data = data([data.group] ~= 2); %%% exclude group
group1 = group(:,1).'; 
% group2 = group(:,2).'; 
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
pain1 = extractfield(data,'Pain1');
pain2 = extractfield(data,'Pain2'); 
%%% Create mask to include subjects with non-missing data
% mask=sum(flag([2,9],:))==2; %%% for Chalder (index1-chldr, index2-Conn)
% Mask for Pain, first line for follow-up pain, second for Conn
mask = ~isnan(pain2); 
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
followUp1 = pain2(mask)'.*group1(mask)'; 
followUp2 = pain2(mask)'.*group3(mask)'; 
% followUp1 = chldr2(mask)'.*group1(mask)'; 
% followUp2 = chldr2(mask)'.*group3(mask)';
covName=['treatment1' 'treatment2' 'gender' 'age' 'site' 'TIV' 'baseline' 'treatment1_followUp' 'treatment2_followUp'];  
% X =[ones(subj-1,1), group1(mask)', group3(mask)']
% X=[group1(mask)',group3(mask)', gender(mask)', age(mask)', site(mask)',TIV(mask)',chldr1(mask)',followUp1,followUp2]; 
X=[group1(mask)',group3(mask)', gender(mask)', age(mask)', site(mask)',TIV(mask)',pain1(mask)',followUp1,followUp2];
% C=[0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1] ; %%% to evaluate whether there is a follow-up fatigue effect in either of the two treatment groups 
C=[0 0 0 0 0 0 0 1 -1]; %%% Check for differences between treatment groups

%%% to create GLM using the covariate and  data
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
    if LME.p_FDR(i) <0.1
        signif{i} = 'YES';
    else
        signif{i} = 'NO'; 
    end
end   
%%% Index for significant connections
Index = find(contains(signif,'YES')); 
