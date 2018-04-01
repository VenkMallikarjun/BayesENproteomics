# BayesENproteomics
Non-linear Bayesian elastic net regression for calculating protein and PTM fold changes from peptide intensities in label-free, bottom-up proteomics on heterogeneous primary human samples. [LINK to PREPRINT WHEN UPLOADED]

# Requirements:
Matlab (2012a or higher) with Bioinformatics and Statistics and Machine Learning toolboxes.

# Usage:
Complete analysis of a dataset using MS1 peptide intensities from a Progenesis QI-formatted .csv spreadsheet (see Progenesis QI folder for examples) can be performed using by calling the BayesENproteomics function. For the mixed species dataset in Fig. 2 and 3 of [PREPRINT], analysis can be performed by calling BayesENproteomics.m as follows:

Human-specific peptide analysis:
>>[HumanProteinOutput,HumanPathwayOutput] = BayesENproteomics('20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
>>'20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv','human',3,donors,{''},'MSCMEF');
>>
Where 'donors' is a row vector = [1,2,3,4,1,2,3,4,1,2,3,4] denoting which donor each MS run is from.

Mouse-specific peptide analysis:
>>[MouseProteinOutput,MousePathwayOutput] = BayesENproteomics('20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
>>'20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv','mouse',3,ones(1,15),{''},'MSCMEF');

Mouse skin technical replicate analysis:
>>TechRepProteinOutput = >>BayesENproteomics('20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv',...
>>'20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv','mouse',3,ones(1,3));


For the PNGase F-treated vs. ctrl samples in Fig 4 of [PREPRINT] can be performed using the peptide list in the Progenesis QI folder (20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv) by calling:

>> PNGaseProteinOutput = BayesENproteomics('20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv',...
>>'20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv','human',2,donors);
>>
Where 'donors' is a row vector = [1,2,3,4,5,1,2,3,4,5] denoting which donor each MS run is from.


