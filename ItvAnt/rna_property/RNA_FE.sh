#!/bin/bash
fn_seq=$(readlink -f $1)
dir_RNA_FE=$2
run_RNA_fold=$3
bn=$(basename $fn_seq .fa)

cd $dir_RNA_FE
mkdir $bn
cd $bn
#predict secondary structure.
$run_RNA_fold -i $fn_seq -o

function get_energ(){ #.fold
    fn=$1
    bn=$(basename $fn .fold)
    nr=$(wc -l $fn|awk '{print $1}')
    case $nr in
        3) E=$(head -n 3 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
            echo -e $bn"\t"$E;;
        6) E1=$(head -n 3 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
            E2=$(head -n 6 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
            Em=$(echo $E1 $E2|awk '{print ($1+$2)/2}')
            echo -e $bn"\t"$Em;;
        9) E1=$(head -n 3 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
           E2=$(head -n 6 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
           E3=$(head -n 9 $fn|tail -n 1|awk '{print $NF}'|sed -e 's/(//g;s/)//g')
           Em=$(echo $E1 $E2|awk '{print ($1+$2+$3)/3}')
            echo -e $bn"\t"$Em;;

    esac
}
for j in *.fold;do
	get_energ $j
    done > ../${bn}.energ
cd ../../
