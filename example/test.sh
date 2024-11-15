#!/bin/bash
gen=small
par=param.txt
mafld=maf.ld.txt

## default simulator
out=./pheno_annot/
mkdir -p ${out}
annot=annot.param.txt
../build/Simulator_annot -g $gen -k 10 -jn 100  -o $out -annot $annot

## simulator with genetic architecture depending on MAF/LD
out=./pheno_mafld/
mkdir -p ${out}
annot=annot.txt
../build/Simulator_mafld -g $gen  -simul_par  $par -maf_ld $mafld  -k 10 -jn 100  -o $out -annot $annot
