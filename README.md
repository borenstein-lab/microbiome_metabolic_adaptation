# THE ASSOCIATIONS BETWEEN PHYSICAL ACTIVITY, MICROBIOME AND METABOLIC ADAPTATION IN SEDENTARY OVERWEIGHT ADULTS

This repository contains the analysis code for the ARE collaboration trial between the Gepner and Borenstein labs. 


## Scripts and analysis notebook overview 

* "analysis": A directory of Rmd notebook generating figures 1-4 and the supplemenrary material. 
* "R": A directory of R scipts used in the analysis:
	1. **plotting_functions.R** - setup of a common plot theme, and templates for box and scatter plots. 
	2. **ml_functions.R** - code for response prediction by univariate logistic regression (including cross validation and result visualization). 
	3. **data_processing_functions.R** - processing of Kraken2 abundance tables
	4. **beta_diversity_functions.R** - generating PCoA plots and running permanova tests. 
