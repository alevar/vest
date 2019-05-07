#!/usr/bin/env python

# ./subtypeRun.py <collection of references> <output directory>
import os
import sys
import subprocess

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

refs={}
writeit=False
def main(args):
	assert os.path.exists(args[1]), "output path does not exist"

	# load references
	curRef=""
	with open(args[0],"r") as f:
		for line in f.readlines():
			if line[0]==">":
				name=line[1:].strip().split(".")[0]
				if name in subtypes:
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

	# save subtypes only
	if not os.path.exists(args[1].rstrip("/")+"/vestDB/"):
		os.makedirs(args[1].rstrip("/")+"/vestDB/")
	outFP=open(args[1].rstrip("/")+"/vestDB/refs.fa","w+")
	for ref in refs:
		outFP.write(">"+ref+"\n")
		outFP.write(refs[ref]+"\n")
	outFP.close()

	# build multiple sequence alignment
	subprocess.call(["muscle",
					 "-in",
					 args[1].rstrip("/")+"/vestDB/refs.fa",
					 "-out",
					 args[1].rstrip("/")+"/vestDB/refs.mus"])

	# build vest database
	subprocess.call(["./vest.py",
					 "build",
					 "-i",
					 args[1].rstrip("/")+"/vestDB/refs.mus",
					 "-o",
					 args[1].rstrip("/")+"/vestDB/msa"])

	#=============================================
	# experiment 1 - simulate sequencing of each RefSeq. Do not remove the subtype only true genome
	#=============================================

	# first create directory and write reads to a file
	if not os.path.exists(args[1].rstrip("/")+"/singles"):
		os.makedirs(args[1].rstrip("/")+"/singles")
	for ref in list(refs):
		outDir=args[1].rstrip("/")+"/"+ref
		if not os.path.exists(outDir):
			os.makedirs(outDir)
			if not os.path.exists(outDir+"/dbs"):
				os.makedirs(outDir+"/dbs")

		outFP=open(outDir+"/bowtie.fa","w+")
		for r in list(refs):
			if not r==ref: # ignore true genome
				outFP.write(">"+r+"\n")
				outFP.write(refs[r]+"\n")
		outFP.close()
		subprocess.call(["bowtie2-build",
						 outDir+"/bowtie.fa",
						 outDir+"/dbs/bowtie"])
		
		outFP=open(outDir+"/true.fa","w+")
		outFP.write(">"+ref+"\n")
		outFP.write(refs[ref]+"\n")
		outFP.close()
		subprocess.call(["bowtie2-build",
						 outDir+"/true.fa",
						 args[1].rstrip("/")+"/singles/"+ref])

		# simulate reads from true genome
		subprocess.call(["wgsim",
						 "-N",
						 "3000",
						 "-r",
						 "0.01",
						 "-1",
						 "151",
						 "-S11",
						 "-d0",
						 "-e0",
						 outDir+"/true.fa",
						 outDir+"/reads.fq","/dev/null"])

	# align to other genomes/subtypes individually
	for ref in list(refs):
		outDir=args[1].rstrip("/")+"/"+ref
		outFP=open(outDir+"/stats.csv","w+")
		outFP.write("al,rate,mult,once,zero\n")

		# align against collection of other subtypes
		res=subprocess.check_output(["bowtie2",
									 "--very-sensitive",
									 "--no-unal",
									 "-p",
									 "4",
									 "-x",
									 outDir+"/dbs/bowtie",
									 "-U",
									 outDir+"/reads.fq",
									 "-S",
									 outDir+"/reads.sam"],stderr=subprocess.STDOUT)
		rate=res.split("\n")[-2].split(" ")[0].rstrip("%")
		mult=res.split("\n")[-3].strip().split(" ")[0]
		once=res.split("\n")[-4].strip().split(" ")[0]
		zero=res.split("\n")[-5].strip().split(" ")[0]
		outFP.write("all,"+rate+","+mult+","+once+","+zero+"\n")

		# realign using Vest
		subprocess.call(["./vest.py",
						 "realign",
						 "-i",
						 outDir+"/reads.sam",
						 "-x",
						 args[1].rstrip("/")+"/vestDB/msa",
						 "-o",
						 outDir+"/reads.realigned"])

		# sort by positions
		subprocess.call(["samtools",
						 "sort",
						 outDir+"/reads.realigned",
						 "-o",
						 outDir+"/reads.realigned.bam"])

		# index for viewing in IGV
		subprocess.call(["samtools",
						 "index",
						 outDir+"/reads.realigned.bam"])

		for r in list(refs):
			if not r==ref: # ignore true genome
				res=subprocess.check_output(["bowtie2",
											 "--very-sensitive",
											 "--no-unal",
											 "-p",
											 "4",
											 "-x",
											 args[1].rstrip("/")+"/singles/"+r,
											 "-U",
											 outDir+"/reads.fq",
											 "-S",
											 outDir+"/"+r+".sam"],stderr=subprocess.STDOUT)

				rate=res.split("\n")[-2].split(" ")[0].rstrip("%")
				mult=res.split("\n")[-3].strip().split(" ")[0]
				once=res.split("\n")[-4].strip().split(" ")[0]
				zero=res.split("\n")[-5].strip().split(" ")[0]
				outFP.write(r+","+rate+","+mult+","+once+","+zero+"\n")

		outFP.close()

if __name__=="__main__":
	main(sys.argv[1:])