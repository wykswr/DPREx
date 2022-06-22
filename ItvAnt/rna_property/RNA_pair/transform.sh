#!/bin/bash
bn=$1
fn_st=./ori_scores/${bn}_pair_seq0.score
fn_ed=./ori_scores/${bn}_pair_seq1.score
awk -F',' '{print $1}' $fn_st|sed 's/CHR/chr/g' > chr
awk -F',' '{print $2}' $fn_st > st
awk -F',' '{print $3}' $fn_st > s1
awk -F',' '{print $2}' $fn_ed > ed
awk -F',' '{print $3}' $fn_ed > s2
paste -d ',' chr st ed s1 s2 > ${bn}_formed.pair
sed -i '1i\chr,st,ed,st_pair,ed_pair' ${bn}_formed.pair
rm chr st ed s1 s2

