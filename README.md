# BayesENproteomics
Non-linear Bayesian elastic net regression for calculating protein and PTM fold changes from peptide intensities in label-free, bottom-up proteomics on heterogeneous primary human samples. 

If you find this method useful please cite the [preprint](https://www.biorxiv.org/content/early/2018/05/10/295527).



## Requirements:
Matlab (2012a or higher) with the Statistics and Machine Learning toolbox. An internet connection is required as data from UniProt and Reactome databases are pulled during analysis.

N.B. All .m and .csv flies to be used in analysis must be placed in your working directory.

## Usage:
Complete analysis of a dataset using MS1 peptide intensities from a Progenesis QI-formatted .csv spreadsheet (see Progenesis QI folder for examples) can be performed using by calling two functions:

### 1. BayesENproteomics.m to perform model fitting and table formatting for proteins and pathways.
### 2. DataProcessEBcontrasts.m for experiments with multiple treatments to perform different comparisions between treatments.


## 1. BayesENproteomics.m
Called as: 

`[ProteinOutput, PathwayOutput] = BayesENproteomics(exp_peps, norm_peps, species, groupnum, donors, ptms, norm_method, mins, pep_fdr, nDB, impute_method);`

Where:
  - ProteinOutput = Structure containing protein- and PTM-level quantification.
  
  - PathwayOutput = Structure containing Reactome pathway-level quantiifcation.
  
  * exp_peps = string containing name and extension (e.g. 'name.csv') of file containing quantification and details for peptides to be used in quantification (organised into Progenesis format; for quantification include only Raw Abundances, do not include Spectral Counts or Normalised Abundances).
  
  * norm_peps = string containing name and extension (e.g. 'name.csv') of file containing quantification and details for peptides to be used for normalisation (organised into Progenesis format). Can be same as exp_peps if normalisation is to be performed against entire dataset. File referenced must have the same number of columns as `exp_peps`.
  
  * species = string denoting species used. Can either be 'mouse', 'human' or a string containing a UniProt proteome ID (e.g. 'UP000002281' for Equus caballus). 
  
  * groupnum = number of experimental treatment conditions, including control.
  
  * donors = numerical row vector containing numbers for each sample denoting which donor a given sample was derived from. E.g. 8 samples from 4 donors, donors = [1,2,3,4,1,2,3,4]; If set to false will assume that each sample comes from a separate unique donor with no pairing between conditions. Enter a vector of ones if donor variability is suspected to be negligible (E.g. for cell lines).
  
  * ptms = Optional. Cell array of strings containing PTMs to look for spelt as they are in exp_peps. E.g. {'Phospho', 'Oxidation'}. Defaults to {''}.
  
  * norm_method = Optional. Character array that determines how peptide intensities are normalised prior to model fitting. 'mean','median' or 'sum' equalises that attribute between all runs. 'MSCMEF' normalises each run in exp_peps to its respective         column in norm_peps. Defaults to 'median'.
  
  * mins = Optional. 2-dimensional vector specifying minimum number of [peptides,proteins] required for a [protein,pathway] to be reported. Defaults to [3,5].
  
  * pep_fdr = Optional. Scalar number >0 and <=1 denotes peptide false discovery p-value cutoff. Following Benjamini-Hochberg adjustment, peptides with identification p-values above this threshold will be discarded. Defaults to 0.2.
  
  * nDB = Optional. Scalar number indicating number of databases used for peptide annotation. Defaults to 1.
  
  * impute_method = Optional. String denoting imputation methods to be used. Options are: 'bpca' - Bayesian PCA decomposition of matrix and reassembly of missing values; 'knn' - mean of 11 nearest neighbours determined by euclidian distance; 'dgd' - random imputation from a downshifted Gaussian distribution; and 'ami' - adaptive multiple imputation as detailed in the preprint. Defaults to 'ami'.


### Mixed-species example

For the mixed species dataset in Fig. 2 and 3 of [Mallikarjun et al. (2018)](https://www.biorxiv.org/content/early/2018/04/06/295527), analysis can be performed by calling BayesENproteomics.m as follows:

Human-specific peptide analysis:

`[HumanProteinOutput,HumanPathwayOutput] = BayesENproteomics...('20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
'20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv','human',3,donors,{''},'MSCMEF');`

Where 'donors' is a row vector = [1,2,3,4,1,2,3,4,1,2,3,4] denoting which donor each MS run is from.


Mouse-specific peptide analysis:

`[MouseProteinOutput,MousePathwayOutput] = BayesENproteomics...('20180319_MouseSkinuniquepeptides_MixedSpecies_BayesENproteomics.csv',...
'20180319_HumanMSCsuniquepeptides_MixedSpecies_BayesENproteomics.csv','mouse',3,ones(1,12),{''},'MSCMEF');`


Mouse skin technical replicate analysis:

`TechRepProteinOutput = BayesENproteomics('20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv',...
'20180319_MouseSkinPeptides_technicalreplicates_BayesENproteomics.csv','mouse',3,ones(1,3));`

### PNGase F-treated example

For the PNGase F-treated vs. ctrl samples in Fig 4 of [Mallikarjun et al. (2018)](https://www.biorxiv.org/content/early/2018/04/06/295527) can be performed using the peptide list in the Progenesis QI folder (20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv) by calling:

`PNGaseProteinOutput = BayesENproteomics('20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv',...
'20180103_MSC_PNGaseFbenchmark_peptidelist_BayesENproteomics.csv','human',2,donors);`

Where 'donors' is a row vector = [1,2,3,4,5,1,2,3,4,5] denoting which donor each MS run is from. 


## 2. DataProcessEBcontrasts.m

For experiments with multiple treatments where different comparisons are necessary, this can be performed calling DataProcessEBcontrasts.m as follows:

`contrasted = DataProcessEBContrasts(ProteinOutput, ctrl, PTM, GroupNum, d0s0);`

Where:
- contrasted = cell array containg ProteinOutput with the specified comparison.

- ProteinOutput = Protein/PTM/Pathway-level output from BayesENproteomics.m (ProteinOutput.Abds/ProteinOutput.PTMs/PathwayOutput.Abds, respectively).

- ctrl = Scalar value. Index of group that all others are to be compared to. E.g. to compare to 2nd group listed in ProteinOutput, ctrl = 2.

- PTM = boolean (true or false). Denotes whether comparison is being performed on a PTM output table (true) or not (false).

- GroupNum = Scalar value. Denotes how many treatments are in the dataset in total.

- d0s0 = Optional nx2 double matrix, where n = GroupNum. Contains shape and scale parameters for chi^2 distrution fit to standard errors for Empirical Bayes comparison. Defaults to value provided in ProteinOutput{1,1}.


For example, for the human-specific protein output produced above, to compare all groups to group 2 rather than group 1 (default):

`HumanProteinOutputAbds_5050ctrl = DataProcessEBContrasts(HumanProteinOutput.Abds, 2, false, 3);`


### Example protein-, PTM- and pathway-level output .csv files used in the creation of figures in [Mallikarjun et al. (2018)](https://www.biorxiv.org/content/early/2018/05/10/295527) can be found in the Analysis Output folder.
