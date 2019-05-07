#!/usr/bin/env bash

hiv="./finalData/HIV1.fa"
hivDB="./finalData/HIV1/HIV1"
msa="./finalData/HIV1.mus"
hxb2="./finalData/HXB2.fa"
hxb2DB="./finalData/HXB2/HXB2"
ann="./finalData/HXB2_features.gff3"
seqRun="/ccb/salz4-2/mpertea/hiv/rna_deep_seq/data/YHo110916"
outDataDir="./out110916"

mkdir -p ./finalData/HIV1
mkdir -p ./finalData/HXB2

# Perform Multiple Sequence Alignment of the HIV1 genomes
# muscle -in ${hiv} -out ${msa} -maxiters 5 -diags

# build HISAT2 index for HXB2 genome
# hisat2-build ${hxb2} ${hxb2DB}

# build HISAT2 index for HIV1 genomes
# hisat2-build ${hiv} ${hivDB}

# perform alignments with hisat and compute statistics
./hisatAlign.sh ${seqRun} ${outDataDir} ${hivDB} ${hxb2DB} 10

# tail -n+2 ./*/stats.csv -q > ./stats.csv