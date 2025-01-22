# Downing-et-al-2020b
Code and data for: Downing PA, Griffin AS, &amp; Cornwallis CK. 2020b. The benefits of help in cooperative birds – non-existent or difficult to detect? American Naturalist, 195: 1085-1091.


Number of supplementary items: three
1. BH_R_Code.R
2. BH_Data_Extraction.txt
2. BH_Tables_S1-S4.xlsx


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: BH_R_Code.R

This R script contains all the code needed to replicate the analyses (including packages and functions).

- Data manipulation (lines 26 to 40)
- Part 1. Helper benefits without phylogeny (lines 64 to 158)
- Part B. Helper benefits with phylogeny (lines 164 to 195)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: BH_Tables_S1-S4.xlsx

This Excel document contains the following sheets

- Table S1 (PRISMA workflow)
	+ systematic search for data on the relationship between helpers and reproductive success
	+ search term: (reproduct* OR group size) AND (binomial species name OR synonym OR English name)
	+ 140 data rows = 140 species
	+ column descriptions:\
		A. common = English name of each species\
		B. species = latin binomial (matches the Jetz et al. nomenclature)\
		C. group type = whether the groups consist of family or non-family members\
		D. n.studies.WofS = number of studies on the species in Web of Science\
		E. additional.studies = other studies on the species (e.g. from the grey literature)\
		F. after.duplicates = number of studies after duplicates removed\
		G. screened = number of studies whose titles / abstract were examined for suitability\
		H. exluded = number of studies excluded based on the first screen\
		I. full.text = number of studies read in full to identify data\
		J. full.text.exluded = number of studies excluded after reading full text\
		K. n.studies = number of studies with suitable data\
		L. n.effect.sizes = number of effect sizes extracted from studies with suitable data

- Table S2:
	+ reproductive success estimates (fledglings) in pairs vs groups for 19 species
	+ read into R (some of the column headings will need changing to match the R code)
	+ column descriptions:\
		A. common = English name of each species\
		B. species = latin binomial (matches the Jetz et al. nomenclature)\
		C. study class = five level factor describing the type of study\
		D. study design = paired vs unpaired\
		E. description = brief overview of the sample sizes and contrasts made\
		F. reproductive success measure = the metric used to measure reproductive success in each study\
		G. source = where the data were extracted from in the study\
		H. r = correlation between helper number and reproductive success\
		I. Zr = Fisher’s Z-transformed correlation coefficient\
		J. Zvar = variance of Zr, calculated as 1/(n-3) where n is the sample size\
		K. n = sample size used to estimate the effect size\
		L. two tailed p value = the p value had the statistical test been two-tailed\
		M. reported p tails = the reported tails in the study\
		N. notes = brief overview of how each effect size was calculated\
		O. breeder quality = how the study controlled for breeder quality\
		P. territory quality = how the study controlled for territory quality\
		Q. reference = study from which the data to calculate effect sizes were extracted

- Table S3
	+ studies which were excluded from the final dataset and the reasons why
	+ 80 data rows
	+ column descriptions:\
		A. common name = English name of each species\
		B. species = latin binomial (matches the Jetz et al. nomenclature)\
		C. reason for exclusion = why the study was not included in the final dataset\
		D. notes = details on the potential sources of data in the study\
		E. reference = the study which was excluded

- Table S4
	+ parameter estimates from statistical models (see BH_R_Code.R for further information)
	+ 36 rows = output from 8 statistical models
	+ variable names:\
		Zr = parameter estimate from the model\
		lwr CI = lower 95% credible interval from the posterior distribution of the model\
		upr CI = upper 95% credible interval from the posterior distribution of the model


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: BH_Data_Extraction.txt

This plain text document contains details of the figures, tables, and text fragments from which data were obtained, and the calculations used to pool data, organised by species.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
