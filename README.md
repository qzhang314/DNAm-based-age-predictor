# DNAm-based-age-predictors
Chronological age predictor based on DNA methylation from Illumina HumanMethylation450/EPIC arrays 

There are two predictors (two sets of coefficients) built based on same set of 13,566 training samples, but using different methods. One is based on Elastic Net, the other is based on best linear unbiased prediction (BLUP). 514 probes have been selected by Elastic Net and BLUP used all (319,607) probes.  

Since our training dataset contains a large number of samples, these samples should not be used as test samples since the prediction results will be biased. Detailed sample information is available in our paper (see Citation).

## How to use the predictor
```console
Rscript pred.R -i input_file -o output_file -a age_file
```

> input_file: a R object file which contains the DNA methylation information for individuals (N * M matrix). N is the number of individuals and M is the number of CpG sites. Beta value is used as DNA methylation measurement.  

> output_file: the name of output file. The output file contains four columns: individual ID, real chronological age, predicted age using predictior based on Elastic Net and predicted age using predictior based on BLUP.

> age_file: an input file which has two column: individual ID and real chronological age. Please note, the first line should be the header.



## Examples 
#### Please use this example to test whether the preditor can work properly in your working environment
```console
Rscript pred.R -i data.rds -o age.pred -a data.age     
```

## Contact us 
If you have any problems about the code and predictor, please contact us: qz2@sanger.ac.uk

## Citation
> Zhang Q, Vallerga C, Walker R, Lin T, Henders A, Montgomery G, He J, Fan D, Fowdar J, Kennedy M, et al: Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing. Genome Medicine, 2019.
