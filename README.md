# BayesENproteomics
Non-linear Bayesian elastic net regression for calculating protein and PTM fold changes from peptide intensities in label-free, bottom-up proteomics on heterogeneous primary human samples. [LINK to PREPRINT WHEN UPLOADED]

## Requirements:
Matlab (2012a or higher) with Bioinformatics and Statistics and Machine Learning toolboxes. An internet connection is required as data from UniProt and Reactome databases is pulled during analysis.

N.B. All .m and .csv flies to be used in analysis must be placed in your working directory.

## Usage:
Complete analysis of a dataset using MS1 peptide intensities from a Progenesis QI-formatted .csv spreadsheet (see Progenesis QI folder for examples) can be performed using by calling two functions:

### 1. BayesENproteomics.m to perform model fitting and table formatting for proteins and pathways.
### 2. DataProcessEBcontrasts.m for experiments with multiple treatments to perform different comparisions between treatments.


### 1. BayesENproteomics.m
Called as: 

`[ProteinOutput, PathwayOutput] = BayesENproteomics(exp_peps, norm_peps, species, groupnum, donors, ptms, norm_method, mins, pep_fdr, nDB);`

Where:
  - ProteinOutput = Structure containing protein- and PTM-level quantification.
  
  - PathwayOutput = Structure containing Reactome pathway-level quantiifcation.
  
  - exp_peps = string containing name and extension (e.g. 'name.csv') of file containing quantification and details for peptides to be used in quantification (organised into Progenesis format).
  
  - norm_peps = string containing name and extension (e.g. 'name.csv') of file containing quantification and details for peptides to be used for normalisation (organised into Progenesis format). Can be same as exp_peps if normalisation is to be performed against entire dataset.
  
  - species = string denoting species used. Currently only 'mouse' or 'human' can be entered.
  
  - groupnum = number of experimental treatment conditions, including control.
  
  - donors = numerical row vector containing numbers for each sample denoting which donor a given sample was derived from. E.g. 8 samples from 4 donors, donors = [1,2,3,4,1,2,3,4]; If set to false will assume that each sample comes from a separate unique donor with no pairing between conditions. Enter a vector of ones if donor variability is suspected to be negligible (E.g. for cell lines or 100% inbred mouse lines).
  
  - ptms = Optional. Cell array of strings containing PTMs to look for spelt as they are in exp_peps. E.g. {'Phospho', 'Oxidation'}.           Defaults to {''}.
  
  - norm_method = Optional. Character array that determines how peptide intensities are normalised prior to model fitting.                     'mean','median' or 'sum' equalises that attribute between all runs. 'MSCMEF' normalises each run in exp_peps to its respective         column in norm_peps. Defaults to 'median'.
  
  - mins = Optional. 2-dimensional vector specifying minimum number of [peptides,proteins] required for a [protein,pathway] to be             reported. Defaults to [3,5].
  
  - pep_fdr = Optional. Scalar number >0 and <=1 denotes peptide false discovery p-value cutoff. Following Benjamini-Hochberg adjustment, peptides with identification p-values above this threshold will be discarded. Defaults to 0.2.
  
  - nDB = Optional. Scalar number indicating number of databases used for peptide annotation. Defaults to 1.
  
For the mixed species dataset in Fig. 2 and 3 of [PREPRINT], analysis can be performed by calling BayesENproteomics.m as follows:

Human-specific peptide analysis:
`[HumanProteinOutput,HumanPathwayOutput] = BayesENproteomics('20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
'20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv','human',3,donors,{''},'MSCMEF');`

Where 'donors' is a row vector = [1,2,3,4,1,2,3,4,1,2,3,4] denoting which donor each MS run is from.


Mouse-specific peptide analysis:

`[MouseProteinOutput,MousePathwayOutput] = BayesENproteomics('20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
'20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv','mouse',3,ones(1,15),{''},'MSCMEF');`


Mouse skin technical replicate analysis:

`TechRepProteinOutput = BayesENproteomics('20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv',...
'20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv','mouse',3,ones(1,3));`


For the PNGase F-treated vs. ctrl samples in Fig 4 of [PREPRINT] can be performed using the peptide list in the Progenesis QI folder (20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv) by calling:

`PNGaseProteinOutput = BayesENproteomics('20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv',...
'20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv','human',2,donors);`

Where 'donors' is a row vector = [1,2,3,4,5,1,2,3,4,5] denoting which donor each MS run is from. Further options for BayesENproteomics can be found in comments within the code.

For experiments with multiple treatments where different comparisons are necessary, this can be performed calling DataProcessEBcontrasts.m. For example, for the human-specific protein output produced above, to compare all groups to group 2 rather than group 1 (default):

`HumanProteinOutputAbds_5050ctrl = DataProcessEBContrasts(HumanProteinOutput.Abds, 2, false, 3);`

Details for DataProcessEBcontrasts.m can be found within the code header. Protein-, PTM- and Pathway-level output .csv files can be found in the Analysis Output folder.
