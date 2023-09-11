gen=small
annot=annot.txt
par=param.txt
mafld=maf.ld.txt
out=./pheno/
../build/Simulator -g $gen  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100    -o $out -annot $annot






