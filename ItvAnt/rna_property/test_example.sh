#!/bin/bash
mkdir fc_test
./main.py --input_file ./example.tsv --output_file fc_test.tsv --dir_rloop_annot /bigdat1/user/fancong/fan_database/Rloop --dir_BER_annot /bigdat1/user/fancong/fan_database/DNA_repair_db --hg19_in /bigdat1/user/fancong/fan_database/hg19.fa --dir_GF /bigdat1/user/fancong/software/GRASP-master --run_Rfold /bigdat1/user/fancong/software/RNAfold/bin/RNAfold --dir_gen_files ./fc_test 
rm fc_test -rf
