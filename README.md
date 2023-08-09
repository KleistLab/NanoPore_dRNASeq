# Sequencing accuracy and systematic errors in nanopore direct RNA sequencing

Code and notebooks to reproduce the results in ONT dRNA-seq accuracy and error analysis

## How to evaluate your own datasets
To compute the general accuracy stats and systematic error profiles from a nanopore dataset, run
```
python eval.py reference.fa aligned.sam
```
where reference.fa is the reference genome/transcriptome in fasta format and aligned.sam is the read alignment from minimap2 (or other aligners). See below for the minimap2 command.

The output will look like this (this example was run on the KO dataset from [Zhong, Z. et al.(2023)](https://doi.org/10.1038/s41467-023-37596-5).
```
Parsing the entire dataset to compute general stats...
The number of sequences in this sample is 910845
100%|█████████████████████████████████████████████████████████████████████████████| 910845/910845 [01:20<00:00, 11325.81it/s]
------------------------------------------------------------
The overall accuracy on this dataset is 0.939 (mismatch 0.018 insertion 0.017 deletion 0.025)
------------------------------------------------------------
Subsampling 200000 reads to compute systematic error profiles...
The number of sequences in this sample is 200000
100%|██████████████████████████████████████████████████████████████████████████████| 200000/200000 [00:48<00:00, 4159.36it/s]
------------------------------------------------------------
The sequencing errors are 52.1% heteropolymers and 47.9% homopolymers
------------------------------------------------------------
The error profiles of 2-mer motifs
         Mis     Ins     Del
GA      0.025   0.014   0.041
AG      0.034   0.009   0.036
GU      0.038   0.019   0.042
UG      0.036   0.013   0.035
GC      0.028   0.012   0.033
CG      0.069   0.018   0.044
AC      0.049   0.013   0.064
CA      0.040   0.015   0.061
AU      0.033   0.021   0.057
UA      0.040   0.014   0.062
UC      0.062   0.015   0.079
CU      0.057   0.024   0.091
------------------------------------------------------------
The confusion matrix of the 4 reference bases ACGU
del     0.036   0.038   0.022   0.042
U       0.011   0.020   0.004   0.930
G       0.007   0.002   0.960   0.002
C       0.004   0.936   0.002   0.017
A       0.942   0.004   0.013   0.009
        A       C       G       U
```

### How to run minimap2 for the correct samfile formatting
To obtain the approriate samfile for the above evaluation, minimap2 should be run with the ```-secondary=no --eqx --sam-hit-only``` options. 
