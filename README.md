## Combining genetics and real-world patient data to fuel ancestry-specific target and drug discovery in Alzheimer’s disease. Yuan Hou, et al, 2023.

A framework illustrating Mendelian Randomization (MR) analysis using Alzheimer’s disease (AD) as a prototypical example. The MR analysis workflow includes two parallel parts. 

# MR-models
The steps of MR-models for drug target identification and validation. The instrumental variables (IVs) corresponding to 1,176 druggable targets were selected from 3 protein quantitative trait loci (pQTL) and 9 expression quantitative trait loci (eQTL) datasets across 5 brain regions. In total, seven GWAS summary statistic datasets with 275,540 AD cases and 1.55 million controls (both African American and European American) were used as outcome datasets. Five MR methods were used to enhance the reproducibility and scientific rigor of drug target findings. 

# Code:
	1_preprocessing.ipynb
	2_MR.r
	3_alzMR_score.ipynb
	4_plot (circos plot in Fig2, Fig3 plot, upset plot..)


# DWAS
A framework illustrating drugome-wide association studies (DWAS). Specifically, four drug cohort designs were systematically evaluated, including exposed vs non exposed model (Exp.), absolute higher exposed vs. low exposed (Abs.), relatively high exposed vs. lowly exposed (Rel.), and medication possession ratio (MPR).

# Code:
	5_DWAS_sample_codes.r
        4_plot (Fig5c plot...)

