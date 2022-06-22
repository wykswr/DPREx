#!/bin/bash
dir_gf_run=$1
dir_RNA_seq=$2
dir_RNA_pair=$3
conda_python3=$4

#prepare,all sequence should be upper-case
for i in ${dir_RNA_seq}/*.fa;do
	bn=$(basename $i .fa)
	cat $i|tr '[a-z]' '[A-Z]' > ${dir_RNA_seq}/${bn}_upper.fa
done

#calculate
test_py=${dir_gf_run}/test.py
mod_fn=${dir_gf_run}/model/PARS_human/ph_25.model
for i in ${dir_RNA_seq}/*_upper.fa;do
	bn=$(basename $i .fa)
	$conda_python3 $test_py --input_fasta=$i --model_file=$mod_fn --window_size=25 --output_file=${dir_RNA_pair}/${bn}.out
done

