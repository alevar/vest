#!/usr/bin/env bash

# samtools view ./data/sorted.bam "JX503073.1:6000-7400" -b > ./data/pJ1.bam && samtools view ./data/sorted.bam "JX503076.1:6000-7400" -b > ./data/pJ2.bam && samtools view ./data/sorted.bam "JX503071.1:6000-7400" -b > ./data/pJ3.bam && samtools view ./data/sorted.bam "JX503080.1:6000-7400" -b > ./data/pJ4.bam && samtools merge ./data/problem.bam ./data/pK.bam ./data/pJ1.bam ./data/pJ2.bam ./data/pJ3.bam ./data/pJ4.bam -f && samtools index ./data/problem.bam

# head -n 11 ./data/problem.sam > ./data/problem2.sam
# samtools view ./data/problem.sam | awk 'NR>=1 && NR <=1029' >> ./data/problem2.sam

./vest.py realign -i ./data/subset.sam -x ./db/subset -o ./res/problem.bam

# mv ./res/problem.bam.sortedName ./res/problem.bam
rm -f ./res/problem.bam.msa.fa.fai && samtools sort ./res/problem.bam > ./res/problem.sorted.bam && samtools index ./res/problem.sorted.bam