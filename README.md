# Biobank-Scale Simulator
Two versions are available:

Simulator_mafld : simulates phenotypes based on a given MAF and LD dependent genetic articheture.

Simulator_annot : simulates phenotypes based on given SNP annotations and their corresponding variance components. 


## Prerequisites
The following packages are required on a linux machine to compile and use the software package.
```
g++
cmake
make
```


```
git clone https://github.com/alipazokit/simulator.git
cd simulator
mkdir build
cd build/
cmake ..
make

```

## Input arguments for  Simulator_mafld 
```
-g :genotype file in bed format
-annot:annotation file in txt format 
  It has M rows (M=number  of SNPs) and single columns. If SNP i is included , then there is  "1" in row i and column j.    Otherwise there is "0"
-o :output directory 
-maf_ld : The file which includes MAF and LD of SNPs. 
-simul_par :
 the  file which contains simulation parameters used in ((Eq.12)  in our paper ) . It has seven columns (percentage of causal SNPs, exponent of LD,exponent of MAF, min(MAF) of causal SNPs, max(MAF) of causal SNPs, total h2, number of simulations).
 -jn : Number of stream blocks. (the more blocks the lower memory usage)
 There is a toy example in the example folder (please take a look at test.sh file ).
 ```
 
## Input arguments for  Simulator_annot
```
-g :genotype file
-annot : annotation file
 
It has M rows (M=number  of SNPs) and K columns (K=number of annotations). If SNP i belongs to annotation j, then there is  "1" in row i and column j.    Otherwise there is "0". (delimiter is " ")

In the first line, please write the true values of variance components.

-o :output directory 
-k : Number of phenotypes
-jn : Number of stream blocks. (the more blocks the lower memory usage)
 
 ```
 
