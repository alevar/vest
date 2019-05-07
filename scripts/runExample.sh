#!/usr/bin/env bash

# run muscle to generate multiple sequence alignment
# muscle -in ./examples/genomes.fa -out ./examples/genomes.mus -maxiters 1 -diags

# build a database and index for the genomes in MSA
./vest.py build -i ./examples/genomes.mus -o ./examples/genomes --draw ./examples/build.png

# realign reads to the MSA
./vest.py realign -i ./examples/reads.sam -x ./examples/genomes -o ./examples/realigned.bam --draw ./examples/realign.png

rm -f ./examples/realigned.bam.msa.fa.fai && samtools sort ./examples/realigned.bam > ./examples/realigned.sorted.bam && samtools index ./examples/realigned.sorted.bam