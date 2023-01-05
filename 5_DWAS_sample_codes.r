#######################################################################################
### Screening analysis

### Define outcome and covariates

AD_idx=data_condition$AD_idx
AD_time=data_condition$AD_time
Cancer=data_condition$Cancer
Cerebrovascular=data_condition$Cerebrovascular
Diabetes=data_condition$Diabetes
Heart=data_condition$Heart
Liver=data_condition$Liver
Pulmonary=data_condition$Pulmonary
Renal=data_condition$Renal
Anxiety=data_condition$Anxiety
Bipolar=data_condition$Bipolar
Depression=data_condition$Depression
Schizophrenia=data_condition$Schizophrenia
Arthritis=data_condition$Arthritis
Fall=data_condition$Fall
HTN=data_condition$HTN
Hearing=data_condition$Hearing
Vision=data_condition$Vision
WgtLoss=data_condition$WgtLoss
SUD=data_condition$Substance
Peripheral=data_condition$Peripheral
Rheumatic=data_condition$Rheumatic
Anemia=data_condition$Anemia
Thyroid=data_condition$Thyroid
Gender=data_condition$Gender
Age=data_condition$Age
Race_B=data_condition$Race_B
Race_H=data_condition$Race_H
Race_A=data_condition$Race_A
Race_U=data_condition$Race_U

### MPR analysis

for(i in 1:ncol(data_drug))
{
 temp_MPR_A=data_drug[,i]
 temp_name=paste("Drug_", i, sep = "")
 assign(temp_name,temp_MPR_A)
}

xnam=paste("Drug_",1:ncol(data_drug), sep="")
fmla=as.formula(paste("Surv(AD_time,AD_idx)~Gender+Age+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+SUD+Peripheral+Rheumatic+Anemia+Thyroid+",paste(xnam, collapse= "+")))
Cox_MPR=coxph(fmla)

### MPR absolute high exposure vs. absolute low exposure analysis

for(i in 1:ncol(data_drug))
{
 temp_MPR=data_drug[,i]
 temp_thres=0.25
 temp_MPR_A=rep(0,length(temp_MPR))
 temp_MPR_A[which(temp_MPR>temp_thres)]=1
 temp_name=paste("Drug_", i, sep = "")
 assign(temp_name,temp_MPR_A)
}

xnam=paste("Drug_",1:ncol(data_drug), sep="")
fmla=as.formula(paste("Surv(AD_time,AD_idx)~Gender+Age+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+SUD+Peripheral+Rheumatic+Anemia+Thyroid+",paste(xnam, collapse= "+")))
Cox_Abs=coxph(fmla)

### MPR relative high exposure vs. relative low analysis

for(i in 1:ncol(data_drug))
{
 temp_MPR=data_drug[,i]
 temp_thres=quantile(temp_MPR[temp_MPR>0],0.25)
 temp_MPR_A=rep(0,length(temp_MPR))
 temp_MPR_A[which(temp_MPR>temp_thres)]=1
 temp_name=paste("Drug_", i, sep = "")
 assign(temp_name,temp_MPR_A)
}

xnam=paste("Drug_",1:ncol(data_drug), sep="")
fmla=as.formula(paste("Surv(AD_time,AD_idx)~Gender+Age+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+SUD+Peripheral+Rheumatic+Anemia+Thyroid+",paste(xnam, collapse= "+")))
Cox_Rel=coxph(fmla)

### MPR exposure vs no-exposure analysis

for(i in 1:ncol(data_drug))
{
 temp_MPR=data_drug[,i]
 temp_thres=0
 temp_MPR_A=rep(0,length(temp_MPR))
 temp_MPR_A[which(temp_MPR>temp_thres)]=1
 temp_name=paste("Drug_", i, sep = "")
 assign(temp_name,temp_MPR_A)
}

xnam=paste("Drug_",1:ncol(data_drug), sep="")
fmla=as.formula(paste("Surv(AD_time,AD_idx)~Gender+Age+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+SUD+Peripheral+Rheumatic+Anemia+Thyroid+",paste(xnam, collapse= "+")))
Cox_Exp=coxph(fmla)

#######################################################################################
### Validation analysis

library(MatchIt)

result_matching_Amlodipine_Exp=matchit(Amlodipine_Exp~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Amlodipine_Exp,method="Nearest",ratio=1)
result_matching_Amlodipine_Abs=matchit(Amlodipine_Abs~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Amlodipine_Abs,method="Nearest",ratio=1)
result_matching_Amlodipine_Rel=matchit(Amlodipine_Rel~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Amlodipine_Rel,method="Nearest",ratio=1)

result_matching_Apixaban_Exp=matchit(Apixaban_Exp~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Apixaban_Exp,method="Nearest",ratio=1)
result_matching_Apixaban_Abs=matchit(Apixaban_Abs~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Apixaban_Abs,method="Nearest",ratio=1)
result_matching_Apixaban_Rel=matchit(Apixaban_Rel~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Apixaban_Rel,method="Nearest",ratio=1)

result_matching_Baclofen_Exp=matchit(Baclofen_Exp~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Baclofen_Exp,method="Nearest",ratio=1)
result_matching_Baclofen_Abs=matchit(Baclofen_Abs~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Baclofen_Abs,method="Nearest",ratio=1)
result_matching_Baclofen_Rel=matchit(Baclofen_Rel~Age+Gender+Race_B+Race_H+Race_A+Race_U+Cerebrovascular+Diabetes+Heart+Liver+Pulmonary+Renal+Anxiety+Bipolar+Depression+Schizophrenia+Arthritis+Fall+HTN+Hearing+Vision+WgtLoss+SUD+Peripheral+Rheumatic+Anemia+Thyroid,data=data_Baclofen_Rel,method="Nearest",ratio=1)

data_matched_Amlodipine_Exp=match.data(result_matching_Amlodipine_Exp)
data_matched_Amlodipine_Abs=match.data(result_matching_Amlodipine_Abs)
data_matched_Amlodipine_Rel=match.data(result_matching_Amlodipine_Rel)

data_matched_Apixaban_Exp=match.data(result_matching_Apixaban_Exp)
data_matched_Apixaban_Abs=match.data(result_matching_Apixaban_Abs)
data_matched_Apixaban_Rel=match.data(result_matching_Apixaban_Rel)

data_matched_Baclofen_Exp=match.data(result_matching_Baclofen_Exp)
data_matched_Baclofen_Abs=match.data(result_matching_Baclofen_Abs)
data_matched_Baclofen_Rel=match.data(result_matching_Baclofen_Rel)

library(survival)

survdiff(Surv(AD_time,AD_idx)~Amlodipine_Exp,data=data_matched_Amlodipine_Exp)
survdiff(Surv(AD_time,AD_idx)~Amlodipine_Abs,data=data_matched_Amlodipine_Abs)
survdiff(Surv(AD_time,AD_idx)~Amlodipine_Rel,data=data_matched_Amlodipine_Rel)

survdiff(Surv(AD_time,AD_idx)~Apixaban_Exp,data=data_matched_Apixaban_Exp)
survdiff(Surv(AD_time,AD_idx)~Apixaban_Abs,data=data_matched_Apixaban_Abs)
survdiff(Surv(AD_time,AD_idx)~Apixaban_Rel,data=data_matched_Apixaban_Rel)

survdiff(Surv(AD_time,AD_idx)~Baclofen_Exp,data=data_matched_Baclofen_Exp)
survdiff(Surv(AD_time,AD_idx)~Baclofen_Abs,data=data_matched_Baclofen_Abs)
survdiff(Surv(AD_time,AD_idx)~Baclofen_Rel,data=data_matched_Baclofen_Rel)

