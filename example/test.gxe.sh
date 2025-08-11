

gen=small
annot=annot.txt
mafld=maf.ld.txt
env=dummy_env.txt
par=param.gxe.txt

## default GxE simulator
out=./pheno_gxe/
mkdir -p ${out}
../build/Simulator_gxe -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100    -o $out -annot $annot

## simulate GxE correlated with E
out=./pheno_gxe_with_corr/
mkdir -p ${out}
../build/Simulator_gxe_with_corr -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100  -o $out -annot $annot -corr_w_E 0.5

## same causal SNPs for G and GxE components, to be saved in causal_snp_file
causal_snp_file=causal_snp_add.txt
out=./pheno_gxe_same_causal/
mkdir -p ${out}
../build/Simulator_gxe_same_causal -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100  -o $out -annot $annot -causal_snp_file ${causal_snp_file}

## specify the additive causal SNPs (G component) in causal_snp_file
causal_snp_file=causal_snp_add.txt
out=./pheno_gxe_causal_from_file/
mkdir -p ${out}
../build/Simulator_gxe_causal_from_file -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100  -o $out -annot $annot -causal_snp_file ${causal_snp_file}

## heavy-tail noise (student's t-dist with dof=4)
out=./pheno_gxe_heavy_tail/
mkdir -p ${out}
../build/Simulator_gxe_heavy_tail -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100  -o $out -annot $annot

## one-hot E
out=./pheno_gxe_one_hot/
par=param.gxe.one_hot.txt
env=dummy_env.bin.txt
mkdir -p ${out}
../build/Simulator_one_hot -g $gen -e $env  -simul_par  $par -maf_ld $mafld   -k 10 -jn 100  -o $out -annot $annot