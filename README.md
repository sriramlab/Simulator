# Biobank-Scale Simulator

## Overview

This repository contains tools for simulating phenotypes at biobank scale, incorporating various genetic and environmental components.

### Simulation Types

1. **Phenotypes with Additive (G) and Noise Components**:
   - **Simulator_annot**: Simulates phenotypes based on SNP annotations and their variance components.
   - **Simulator_mafld**: Simulates phenotypes based on MAF and LD-dependent genetic architecture.

2. **Phenotypes with Additive (G), Gene-by-Environment (GxE), Noise-by-Environment (NxE), and Noise Components**:
   - **Simulator_gxe**: Simulates phenotypes with G + GxE + NxE components (or combinations as specified in the parameter file).
   - **Simulator_gxe_with_corr**: Simulates phenotypes correlated with the environment.
   - **Simulator_gxe_same_causal**: Simulates phenotypes with the same causal SNPs for G and GxE components.
   - **Simulator_gxe_causal_from_file**: Simulates phenotypes with user-specified causal SNPs for the G component.
   - **Simulator_gxe_heavy_tail**: Simulates phenotypes with noise from a heavy-tailed distribution (Student's t-distribution with DoF = 4).

3. **Phenotypes with Additive (G), Gene-by-gene (GxG), and Noise Components**: 
    - **Simulator_gxg**: Simulates phenotypes with G + GxG components.

---

## Prerequisites

The following packages are required on a Linux machine to compile and use this software:

```
g++
cmake
make
```

## Installation

To set up the simulator, clone the repository and build the project:

```bash
git clone https://github.com/sriramlab/Simulator.git
cd Simulator
mkdir build
cd build/
cmake ..
make
```

## Inputs

### Simulator_annot
#### Parameters:

* `-g`: Genotype file in BED format
* `-annot`: Annotation file with M+1 rows and K columns (M = number of SNPs, K = number of annotations).
  * SNP annotations are binary (1 = belongs to annotation, 0 = otherwise).
  * The first line specifies the true variance components.
* `-o`: Output directory
* `-k`: Number of phenotypes
* `-jn`: Number of stream blocks (to reduce memory usage).

### Simulator_mafld
#### Parameters:

* `-g`: Genotype file 
* `-annot`: Annotation file
  * M rows (M = number of SNPs), 1 column (1 = included, 0 = not included).
* `-o`: Output directory
* `-maf_ld`: File with MAF and LD of SNPs.
* `-simul_par`: Simulation parameters file with 7 columns.
  * Seven parameters: percentage of causal SNPs, exponent of LD,exponent of MAF, min(MAF) of causal SNPs, max(MAF) of causal SNPs, total h2, number of simulations.
  * Refer to this [publication](https://doi.org/10.1038/s41467-020-17576-9) for details.
* `-jn`: Number of stream blocks (to reduce memory usage).


### Simulator_gxe
#### Parameters:

* `-g`: Genotype file in BED format
* `-e`: Environment file in TXT format
  * N rows (N = number of individuals), L columns (L = number of environments).
* `-annot`: Annotation file in TXT format
* `-o`: Output directory
* `-maf_ld`: File with MAF and LD of SNPs
* `-simul_par`: Simulator parameter file with 3 lines (for G, GxE, and NxE components).
* `-jn`: Number of stream blocks

### Simulator_gxe_with_corr
#### Parameters: 
* Same as **Simulator_gxe** plus:
  * `-corr_w_E`: Correlation with the environment (E).

### Simulator_gxe_same_causal
#### Parameters: 
* Same as **Simulator_gxe** plus:
  * `-causal_snp_file`: Output file specifying the causal SNPs (1 = causal, 0 = not causal).

### Simulator_gxe_same_causal
#### Parameters: 
* Same as **Simulator_gxe** plus:
  * `-causal_snp_file`: Input file specifying causal SNPs for the additive (G) component.


### Simulator_gxe_heavy_tail
#### Parameters: 
* Same as **Simulator_gxe** 

### Simulator_gxg
#### Parameters:
* `-g` : Genotype file in BED format
* `-annot`: Annotation file in TXT format
  * M rows (M=number  of SNPs) and 3 columns.
  * First column corresponds to SNPs with additive effect. 
  * Second and third column correspond to two groups of SNPs which have interaction with each other (there is no interaction inside each group). Third column/group is assumed to be smaller than second.
  * If SNP i belongs to annotation j, then there is  "1" in row i and column j. Otherwise there is "0". (delimiter is " ")
* `-jn`: The number of stream blocks. (the higher number of blocks, the lower memory usage)
* `-o`: The path of output file.


## Examples

Toy examples are provided in the `example` folder:
* `test.sh`
* `test.gxe.sh`
* `test.gxg.sh`


## Citations

1. Ali Pazokitoroudi, Zhengtong Liu, Andrew Dahl, Noah Zaitlen, Saharon Rosset, Sriram Sankararaman.
   *AJHG (2024)*; DOI: [10.1016/j.ajhg.2024.05.015](https://doi.org/10.1016/j.ajhg.2024.05.015)

2. Ali Pazokitoroudi, Yue Wu, Kathryn S. Burch, Kangcheng Hou, Aaron Zhou, Bogdan Pasaniuc, Sriram Sankararaman.
   *Nature Communications (2020)*; DOI: [10.1038/s41467-020-17576-9](https://doi.org/10.1038/s41467-020-17576-9)

3. Boyang Fu*, Ali Pazokitoroudi*, Albert Xue, Aakarsh Anand, Prateek Anand, Noah Zaitlen, Sriram Sankararaman.
   *bioRxiv*, 2020. DOI: [10.1101/2023.09.10.557084](https://doi.org/10.1101/2023.09.10.557084)
