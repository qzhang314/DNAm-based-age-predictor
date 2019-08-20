# DNAm-based-age-predictors
Chronological age predictor based on DNA methylation from Illumina HumanMethylation450/EPIC arrays 

There are two predictors (two sets of coefficients) built based on same set of 13,566 training samples, but using different methods. One is based on Elastic Net, the other is based on best linear unbiased prediction (BLUP). 514 probes have been selected by Elastic Net and BLUP used all (319,607) probes.  


## how to use the predictor
Rscript pred.R -i input_file -o output_file -a age_file

> input_file: a R object file which contains the DNA methylation information for individuals (N * M matrix). N is the number of individuals and M is the number of CpG sites. Beta value is used as DNA methylation measurement.  

> output_file: the name of output file. The output file contains four columns: individual ID, real chronological age, predicted age using predictior based on Elastic Net and predicted age using predictior based on BLUP.

> age_file: an input file which has two column: individual ID and real chronological age. Please note, the first line should be the header.

## Examples 
###### please use this example to test whether the preditor can work properly in your working environment
Rscript pred.R data.rds age.pred data.age     


## Contact us 
If you have any problems about the code and predictor, please contact us: q.zhang@uq.edu.au

## Citation
