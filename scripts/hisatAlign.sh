#!/usr/bin/env bash

# ./hisatAlign.sh /ccb/salz4-2/mpertea/hiv/rna_deep_seq/data/YHo102816 ./out102816 ./finalData/refs/HIV1/HIV1 ./finalData/HXB2/HXB2 10 /ccb/salz3-scratch/avaraby/TrimGalore-0.4.3/trim_galore

inputDir=${1}
outputDir=${2}
hivDB=${3}
hxb2DB=${4}
threads=${5}

mkdir -p ${outputDir}

touch ${outputDir}/stats.csv
echo run,sample,R,hiv,hxb2 > ${outputDir}/stats.csv

for file in ${inputDir}/*R1*fastq.gz ; do
	sampleR1Base=$(basename ${file})
	sampleR1="${sampleR1Base%.*.*}"
	sample="${sampleR1%_R1*}"
	sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz
	baseEnd="${sampleR1##*_R1}"

	TOTAL_TIME=0

	echo "========================================"
	echo ${sample}
	echo "========================================"

	echo ADAPTOR AND QUALITY TRIMMING R1
	trim_galore -q 5 \
					--phred33 \
					--length 45 \
					-o ./${outputDir}/ \
					--dont_gzip ${inputDir}/${sampleR1Base}
	trimmedFQ=./${outputDir}/"${sampleR1%_R1*}"_R1"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}/${sample}${baseEnd}.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R1 TO ALL HIV1
	hisat2 -p ${threads} \
			--very-sensitive \
			--no-unal \
			--phred33 \
			-x ${hivDB} \
			-U ${outputDir}/${sample}${baseEnd}.fastq \
			-S ${outputDir}/${sample}${baseEnd}_R1.hiv.sam \
			--al ${outputDir}/${sample}${baseEnd}.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R1 TO HXB2
	hisat2 -p ${threads} \
			--very-sensitive \
			--no-unal \
			--phred33 \
			-x ${hxb2DB} \
			-U ${outputDir}/${sample}${baseEnd}.al.fastq \
			-S ${outputDir}/${sample}${baseEnd}_R1.hxb2.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	rm -f ${outputDir}/${sample}${baseEnd}.al.fastq
	rm -f ${outputDir}/${sample}${baseEnd}.fastq

	hiv="$(samtools view ${outputDir}/${sample}${baseEnd}_R1.hiv.sam | cut -f1 -d$'\t' | sort | uniq | wc -l)"
	hxb2="$(samtools view ${outputDir}/${sample}${baseEnd}_R1.hxb2.sam | cut -f1 -d$'\t' | sort | uniq | wc -l)"
	echo  ${outputDir},${sample}${baseEnd},R1,${hiv},${hxb2} >> ${outputDir}/stats.csv

	#===============R2================

	echo ADAPTOR AND QUALITY TRIMMING R2
	trim_galore -q 5\
					-phred33 \
					--length 45 \
					-o ./${outputDir}/ \
					--dont_gzip ${inputDir}/${sampleR2Base}
	trimmedFQ=./${outputDir}/"${sampleR1%_R1*}"_R2"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}/${sample}${baseEnd}.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R2 TO ALL HIV1
	hisat2 -p ${threads} \
			--very-sensitive \
			--no-unal \
			--phred33 \
			-x ${hivDB} \
			-U ${outputDir}/${sample}${baseEnd}.fastq \
			-S ${outputDir}/${sample}${baseEnd}_R2.hiv.sam \
			--al ${outputDir}/${sample}${baseEnd}.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R2 TO HXB2
	hisat2 -p ${threads} \
			--very-sensitive \
			--no-unal \
			--phred33 \
			-x ${hxb2DB} \
			-U ${outputDir}/${sample}${baseEnd}.al.fastq \
			-S ${outputDir}/${sample}${baseEnd}_R2.hxb2.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	rm -f ${outputDir}/${sample}${baseEnd}.al.fastq
	rm -f ${outputDir}/${sample}${baseEnd}.fastq


	T_DUR="$(($TOTAL_TIME / 60)) minutes and $(($TOTAL_TIME % 60)) seconds"
	echo TOTAL TIME ELAPSED: ${T_DUR}
	echo "$(date)"

	# now that this is done, output the statistic into a csv
	hiv="$(samtools view ${outputDir}/${sample}${baseEnd}_R2.hiv.sam | cut -f1 -d$'\t' | sort | uniq | wc -l)"
	hxb2="$(samtools view ${outputDir}/${sample}${baseEnd}_R2.hxb2.sam | cut -f1 -d$'\t' | sort | uniq | wc -l)"
	echo  ${outputDir},${sample}${baseEnd},R2,${hiv},${hxb2} >> ${outputDir}/stats.csv
done