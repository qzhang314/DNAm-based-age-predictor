# DNAm-based-age-predictor
Chronological age predictor based on DNA methylation from Illumina HumanMethylation450/EPIC arrays 

There are two predictors (two sets of coefficients) built based on same set of 13,566 training samples, but using different methods. One is based on Elastic Net, the other is based on best linear unbiased prediction (BLUP). 

################## how to use the predictor ###############

Rscript pred.R input_file output_file age_file


################## one example ############################

Rscript pred.R data.rds age.pred data.age     #### please use this example to test whether the preditor can work properly in your working environment


Contact us: If you have any problems about the code and predictor, please contact us: q.zhang@uq.edu.au
