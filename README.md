# Normal tissue complication probability (NTCP) modelling using functional data analysis methods

Code for radiotherapy NTCP modelling using functional data analysis methods.

To perform the modelling use the ntcpModellingFunctionalDataAnalysis.R R script.

Values for the clinical covariates and toxicity outcomes should be entered into the [toxicityName]FDA.csv (e.g. mucositisFDA.csv) file. The dose metric histogram data should be enetered into the [organ-at-risk][doseMetric]ForFDA.csv (e.g. OMdvhForFDA.csv) file.

The following dependencies should be installed:

- R
- The R packages fda.usc, ggplot2, caret, reshape2, rms, corrplot, CalibrationCurves, glmnet

Details of the modelling can be found in the article:

Functional data analysis applied to modeling of severe acute mucositis and dysphagia resulting from head and neck radiation therapy. Dean  JA et al. Int J Radiat Oncol Biol Phys 2016.

Please consider citing the above article in publications using the code.

This script is intended for research purposes only.

If you encounter any problems trying to use the code or wish to discuss any aspect of it please feel free to contact me at jamie.adam.dean@gmail.com.

