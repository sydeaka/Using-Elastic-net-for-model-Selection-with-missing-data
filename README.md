# Using-Elastic-net-for-model-Selection-with-missing-data

In this script, I demonstrate the use of elastic net to facilitate
model selection given a large number of highly correlated candidate covariates
with some covariate data missing.

Details about the method are provided in published manuscript, 
 "Survival in pulmonary arterial hypertension patients awaiting lung transplantation.",
  http://www.ncbi.nlm.nih.gov/pubmed/24074527 

Written for R version 3.2.1, but will likely work for other versions.

Note: Given a high performance computing platform, it is trivial to 
 run the 12 computations (4 imputations x 3 values of alpha) in parallel.
 The example is written so that it can be easily run on a single computer with
 4 processors, even on a Windows machine.

Author: Sydeaka Watson, PhD
Institution: The University of Chicago, Biostatistics Laboratory
Contact: sydeakawatson@gmail.com
