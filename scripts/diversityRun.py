#!/usr/bin/env python

# this script examines the effect on alignment rate and accuracy of true genome reconstruction
# with an increase in the number of sequences present in the alignment

import os
import sys
import subprocess
import random

sA1=["U51190","AF004885","AF069670","AF484509"] # not present: AF484509, AF069670
sA2=["AF286238","AF286237"]
sB=["K03455","AY173951","AY423387","AY331295"] # not present: AY331295
sC=["U46016","U52953","AF067155","AY772699"]
sD=["K03454","U88824","AY371157","AY253311"] # not present: AY371157, AY253311
sF1=["AF005494","AF077336","AF075703","AJ249238"] # not present: AJ249238
sF2=["AJ249236","AJ249237","AY371158","AF377956"] # not present: AY371158, AF377956, AJ249237, AJ249236
sG=["AF061642","AF061640","AF084936","U88826"]
sH=["AF005496","AF190127","AF190128"] # not present: AF190128
sJ=["AF082394","AF082395"]
sK=["AJ249235","AJ249239"] # not present: AJ249235, AJ249239

subtypes=sA1+sA2+sB+sC+sD+sF1+sF2+sG+sH+sJ+sK
subRefs={}

refs={}

def main(args):
	assert os.path.exists(args[1]), "output path does not exist"
	assert os.path.exists(args[2]), "reads file does not exist"
	assert os.path.exists(args[3]+'.g.npy'), "graph file does not exist"
	assert os.path.exists(args[3]+'.v.npy'), "vertex file does not exist"
	assert os.path.exists(args[3]+'.db'), "vertex file does not exist"

	# begin by building a database of just the available subtype genomes
	curRef=""
	with open(args[0],"r") as f:
		for line in f.readlines():
			if line[0]==">":
				name=line[1:].strip().split(".")[0]
				if name in subtypes:
					subRefs[name]=""
					curRef=name
					writeit=True
				else:
					writeit=False
			elif not curRef=="" and writeit:
				subRefs[curRef]=subRefs[curRef]+line.strip()
			elif not writeit:
				continue
			else:
				print("something wrong with the reference")

	# # save subtypes only as a single database
	# if not os.path.exists(args[1].rstrip("/")+"/vestDB/"):
	# 	os.makedirs(args[1].rstrip("/")+"/vestDB/")
	# outFP=open(args[1].rstrip("/")+"/vestDB/refs.fa","w+")
	# for sub in subRefs:
	# 	outFP.write(">"+sub+"\n")
	# 	outFP.write(subRefs[sub]+"\n")
	# outFP.close()

	# # build multiple sequence alignment
	# subprocess.call(["muscle",
	# 				 "-in",
	# 				 args[1].rstrip("/")+"/vestDB/refs.fa",
	# 				 "-out",
	# 				 args[1].rstrip("/")+"/vestDB/refs.mus"])

	# # build vest database
	# subprocess.call(["./vest.py",
	# 				 "build",
	# 				 "-i",
	# 				 args[1].rstrip("/")+"/vestDB/refs.mus",
	# 				 "-o",
	# 				 args[1].rstrip("/")+"/vestDB/msa"])

	# # build a bowtie database for the collection of subtypes
	# subprocess.call(["bowtie2-build",
	# 				 args[1].rstrip("/")+"/vestDB/refs.fa",
	# 				 args[1].rstrip("/")+"/vestDB/refs"])

	# outFP=open(args[1].rstrip("/")+"/stats.csv","w+")
	# res=subprocess.check_output(["bowtie2",
	# 							 "--very-sensitive",
	# 							 "--no-unal",
	# 							 "-p",
	# 							 "4",
	# 							 "-x",
	# 							 args[1].rstrip("/")+"/vestDB/refs",
	# 							 "-U",
	# 							 args[2],
	# 							 "-S",
	# 							 args[1].rstrip("/")+"/"+"/reads_base.sam"],stderr=subprocess.STDOUT)

	# rate=res.split("\n")[-2].split(" ")[0].rstrip("%")
	# mult=res.split("\n")[-3].strip().split(" ")[0]
	# once=res.split("\n")[-4].strip().split(" ")[0]
	# zero=res.split("\n")[-5].strip().split(" ")[0]
	# outFP.write("0,"+rate+","+mult+","+once+","+zero+"\n")

	# # realign using Vest
	# subprocess.call(["./vest.py",
	# 				 "realign",
	# 				 "-i",
	# 				 args[1].rstrip("/")+"/"+"/reads_base.sam",
	# 				 "-x",
	# 				 args[1].rstrip("/")+"/vestDB/msa",
	# 				 "-o",
	# 				 args[1].rstrip("/")+"/"+"/reads.realigned_base"])

	# # sort by positions
	# subprocess.call(["samtools",
	# 				 "sort",
	# 				 args[1].rstrip("/")+"/"+"/reads.realigned_base",
	# 				 "-o",
	# 				 args[1].rstrip("/")+"/"+"/reads.realigned_base.bam"])

	# # index for viewing in IGV
	# subprocess.call(["samtools",
	# 				 "index",
	# 				 args[1].rstrip("/")+"/"+"/reads.realigned_base.bam"])

	# next build a collection of complete genomes excluding those in subtypes
	curRef=""
	with open(args[0],"r") as f:
		for line in f.readlines():
			if line[0]==">":
				name=line[1:].strip().split(".")[0]
				if not name in subtypes and "complete genome" in line and not "unverified" in line.lower() and not "RNA" in line  and not "CRF" in line:
					refs[name]=""
					curRef=name
					writeit=True
				else:
					writeit=False
			elif not curRef=="" and writeit:
				refs[curRef]=refs[curRef]+line.strip()
			elif not writeit:
				continue
			else:
				print("something wrong with the reference")

	# now make sure those are complete genomes:
	for ref in list(refs):
		if not len(set(refs[ref]).difference(set("ACTG")))==0 and len(refs[ref])>9000 and len(refs[ref])<10000:
			del refs[ref]

	# # next subsample those based on the accession number (AA###)
	# newRefs={}
	# tmpNames={}
	# curAcc=""
	# for ref in sorted(list(refs)):
	# 	accession=ref[:5]
	# 	if not accession in tmpNames:
	# 		if not curAcc=="": # last acession is over, now randomly pick one reference from the collection of related
	# 			pick=random.choice(tmpNames[curAcc])
	# 			newRefs[pick]=refs[pick]
	# 		tmpNames[accession]=[ref]
	# 		curAcc=accession
	# 	else:
	# 		tmpNames[accession].append(ref)

	# newRefsList=list(newRefs)
	# random.shuffle(newRefsList)
	# for i in range(0,len(newRefsList),20):
	# 	for ref in newRefsList[i:i+20]:
	# 		subRefs[ref]=refs[ref]

	# 	# build database, msa, alignment, etc
	# 	outFP2=open(args[1].rstrip("/")+"/vestDB/refs_"+str(i)+".fa","w+")
	# 	for sub in subRefs:
	# 		outFP2.write(">"+sub+"\n")
	# 		outFP2.write(subRefs[sub]+"\n")
	# 	outFP2.close()

	# 	# build multiple sequence alignment
	# 	subprocess.call(["muscle",
	# 					 "-in",
	# 					 args[1].rstrip("/")+"/vestDB/refs_"+str(i)+".fa",
	# 					 "-out",
	# 					 args[1].rstrip("/")+"/vestDB/refs_"+str(i)+".mus"])

	# 	# build vest database
	# 	subprocess.call(["./vest.py",
	# 					 "build",
	# 					 "-i",
	# 					 args[1].rstrip("/")+"/vestDB/refs_"+str(i)+".mus",
	# 					 "-o",
	# 					 args[1].rstrip("/")+"/vestDB/msa_"+str(i)])

	# 	# build a bowtie database for the collection of subtypes
	# 	subprocess.call(["bowtie2-build",
	# 					 args[1].rstrip("/")+"/vestDB/refs_"+str(i)+".fa",
	# 					 args[1].rstrip("/")+"/vestDB/refs_"+str(i)])

	# 	res=subprocess.check_output(["bowtie2",
	# 								 "--very-sensitive",
	# 								 "--no-unal",
	# 								 "-p",
	# 								 "4",
	# 								 "-x",
	# 								 args[1].rstrip("/")+"/vestDB/refs_"+str(i),
	# 								 "-U",
	# 								 args[2],
	# 								 "-S",
	# 								 args[1].rstrip("/")+"/"+"/reads_"+str(i)+".sam"],stderr=subprocess.STDOUT)

	# 	rate=res.split("\n")[-2].split(" ")[0].rstrip("%")
	# 	mult=res.split("\n")[-3].strip().split(" ")[0]
	# 	once=res.split("\n")[-4].strip().split(" ")[0]
	# 	zero=res.split("\n")[-5].strip().split(" ")[0]
	# 	outFP.write(str(i)+","+rate+","+mult+","+once+","+zero+"\n")

	# 	# realign using Vest
	# 	subprocess.call(["./vest.py",
	# 					 "realign",
	# 					 "-i",
	# 					 args[1].rstrip("/")+"/"+"/reads_"+str(i)+".sam",
	# 					 "-x",
	# 					 args[1].rstrip("/")+"/vestDB/msa_"+str(i),
	# 					 "-o",
	# 					 args[1].rstrip("/")+"/"+"/reads.realigned_"+str(i)])

	# 	# sort by positions
	# 	subprocess.call(["samtools",
	# 					 "sort",
	# 					 args[1].rstrip("/")+"/"+"/reads.realigned_"+str(i),
	# 					 "-o",
	# 					 args[1].rstrip("/")+"/"+"/reads.realigned_"+str(i)+".bam"])

	# 	# index for viewing in IGV
	# 	subprocess.call(["samtools",
	# 					 "index",
	# 					 args[1].rstrip("/")+"/"+"/reads.realigned_"+str(i)+".bam"])

	# outFP.close()

	#================================================================
	# COMPLETE DATABASE
	#================================================================
	# lastly build a complete genome bowtie database and align to that
	refs.update(subRefs)
	outFP2=open(args[1].rstrip("/")+"/vestDB/refs_complete.fa","w+")
	with open(args[0],"r") as hivFP:
		for line in hivFP.readlines():
			if line[0]==">":
				fullName=line[1:].strip()
				name=line[1:].strip().split(".")[0]
				if name in refs:
					outFP2.write(">"+fullName+"\n")
					outFP2.write(refs[name]+"\n")
	outFP2.close()

	# build a bowtie database for the collection of subtypes
	subprocess.call(["bowtie2-build",
					 args[1].rstrip("/")+"/vestDB/refs_complete.fa",
					 args[1].rstrip("/")+"/vestDB/refs_complete"])

	res=subprocess.check_output(["bowtie2",
								 "--very-sensitive",
								 "--no-unal",
								 "-p",
								 "4",
								 "-x",
								 args[1].rstrip("/")+"/vestDB/refs_complete",
								 "-U",
								 args[2],
								 "-S",
								 args[1].rstrip("/")+"/"+"/reads_complete.sam"],stderr=subprocess.STDOUT)

	rate=res.split("\n")[-2].split(" ")[0].rstrip("%")
	mult=res.split("\n")[-3].strip().split(" ")[0]
	once=res.split("\n")[-4].strip().split(" ")[0]
	zero=res.split("\n")[-5].strip().split(" ")[0]
	print("complete,"+rate+","+mult+","+once+","+zero+"\n")

	# realign using Vest
	subprocess.call(["./vest.py",
					 "realign",
					 "-i",
					 args[1].rstrip("/")+"/"+"/reads_complete.sam",
					 "-x",
					 args[3],
					 "-o",
					 args[1].rstrip("/")+"/"+"/reads.realigned_complete"])

	# sort by positions
	subprocess.call(["samtools",
					 "sort",
					 args[1].rstrip("/")+"/"+"/reads.realigned_complete",
					 "-o",
					 args[1].rstrip("/")+"/"+"/reads.realigned_complete.bam"])

	# index for viewing in IGV
	subprocess.call(["samtools",
					 "index",
					 args[1].rstrip("/")+"/"+"/reads.realigned_complete.bam"])



if __name__=="__main__":
	main(sys.argv[1:])