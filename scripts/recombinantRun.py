#!/usr/bin/env python

# ./subtypeRun.py <collection of references> <output directory>
import os
import sys
import subprocess

# recombinants are marked with CRF in the name

# this script extracts all recombinants
# aligns them to each subtype
# aligns them against all subtypes together
# compares the results of these alignments

sA1 = ["U51190", "AF004885", "AF069670", "AF484509"]  # not present: AF484509, AF069670
sA2 = ["AF286238", "AF286237"]
sB = ["K03455", "AY173951", "AY423387", "AY331295"]  # not present: AY331295
sC = ["U46016", "U52953", "AF067155", "AY772699"]
sD = ["K03454", "U88824", "AY371157", "AY253311"]  # not present: AY371157, AY253311
sF1 = ["AF005494", "AF077336", "AF075703", "AJ249238"]  # not present: AJ249238
sF2 = ["AJ249236", "AJ249237", "AY371158", "AF377956"]  # not present: AY371158, AF377956, AJ249237, AJ249236
sG = ["AF061642", "AF061640", "AF084936", "U88826"]
sH = ["AF005496", "AF190127", "AF190128"]  # not present: AF190128
sJ = ["AF082394", "AF082395"]
sK = ["AJ249235", "AJ249239"]  # not present: AJ249235, AJ249239

subtypes = sA1 + sA2 + sB + sC + sD + sF1 + sF2 + sG + sH + sJ + sK
subRefs = {}

refs = {}
writeit = False


def main(args):
    assert os.path.exists(args[1]), "output path does not exist"

    # load references and select recombinant genomes
    curRef = ""
    with open(args[0], "r") as f:
        for line in f.readlines():
            if line[0] == ">":
                if " CRF" in line[1:].strip() and "complete" in line[1:].strip():
                    name = line[1:].strip().split(".")[0]
                    refs[name] = ""
                    curRef = name
                    writeit = True
                else:
                    writeit = False
            elif not curRef == "" and writeit:
                refs[curRef] = refs[curRef] + line.strip()
            elif not writeit:
                continue
            else:
                print("something wrong with the reference")

    # now make sure those are complete genomes:
    for ref in list(refs):
        if not len(set(refs[ref]).difference(set("ACTG"))) == 0:
            del refs[ref]

    # now load, extract and save subtypes
    curRef = ""
    with open(args[0], "r") as f:
        for line in f.readlines():
            if line[0] == ">":
                name = line[1:].strip().split(".")[0]
                if name in subtypes:
                    subRefs[name] = ""
                    curRef = name
                    writeit = True
                else:
                    writeit = False
            elif not curRef == "" and writeit:
                subRefs[curRef] = subRefs[curRef] + line.strip()
            elif not writeit:
                continue
            else:
                print("something wrong with the reference")

    # save subtypes only as a single database
    if not os.path.exists(args[1].rstrip("/") + "/vestDB/"):
        os.makedirs(args[1].rstrip("/") + "/vestDB/")
    outFP = open(args[1].rstrip("/") + "/vestDB/refs.fa", "w+")
    for sub in subRefs:
        outFP.write(">" + sub + "\n")
        outFP.write(subRefs[sub] + "\n")
    outFP.close()

    # build multiple sequence alignment
    subprocess.call(["muscle",
                     "-in",
                     args[1].rstrip("/") + "/vestDB/refs.fa",
                     "-out",
                     args[1].rstrip("/") + "/vestDB/refs.mus"])

    # build vest database
    subprocess.call(["./vest.py",
                     "build",
                     "-i",
                     args[1].rstrip("/") + "/vestDB/refs.mus",
                     "-o",
                     args[1].rstrip("/") + "/vestDB/msa"])

    # build a bowtie database for the collection of subtypes
    subprocess.call(["bowtie2-build",
                     args[1].rstrip("/") + "/vestDB/refs.fa",
                     args[1].rstrip("/") + "/vestDB/refs"])

    # now create individual subtype databases
    if not os.path.exists(args[1].rstrip("/") + "/singles"):
        os.makedirs(args[1].rstrip("/") + "/singles")
    for sub in subRefs:
        outFP = open(args[1].rstrip("/") + "/singles/" + sub + ".fa", "w+")
        outFP.write(">" + sub + "\n")
        outFP.write(subRefs[sub] + "\n")
        outFP.close()
        subprocess.call(["bowtie2-build",
                         args[1].rstrip("/") + "/singles/" + sub + ".fa",
                         args[1].rstrip("/") + "/singles/" + sub])

    # ================================
    # now onto the actual experiment
    # ================================

    # first create directory and write reads to a file
    for ref in list(refs):
        outDir = args[1].rstrip("/") + "/" + ref
        if not os.path.exists(outDir):
            os.makedirs(outDir)
            if not os.path.exists(outDir + "/dbs"):
                os.makedirs(outDir + "/dbs")

        outFP = open(outDir + "/true.fa", "w+")
        outFP.write(">" + ref + "\n")
        outFP.write(refs[ref] + "\n")
        outFP.close()

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
                         outDir + "/true.fa",
                         outDir + "/reads.fq", "/dev/null"])

        # now align to the collection of subtypes
        res = subprocess.check_output(["bowtie2",
                                       "--very-sensitive",
                                       "--no-unal",
                                       "-p",
                                       "4",
                                       "-x",
                                       args[1].rstrip("/") + "/vestDB/refs",
                                       "-U",
                                       outDir + "/reads.fq",
                                       "-S",
                                       outDir + "/reads.sam"], stderr=subprocess.STDOUT)

        outFP = open(outDir + "/stats.csv", "w+")
        rate = res.split("\n")[-2].split(" ")[0].rstrip("%")
        mult = res.split("\n")[-3].strip().split(" ")[0]
        once = res.split("\n")[-4].strip().split(" ")[0]
        zero = res.split("\n")[-5].strip().split(" ")[0]
        outFP.write("all," + rate + "," + mult + "," + once + "," + zero + "\n")

        # realign using Vest
        subprocess.call(["./vest.py",
                         "realign",
                         "-i",
                         outDir + "/reads.sam",
                         "-x",
                         args[1].rstrip("/") + "/vestDB/msa",
                         "-o",
                         outDir + "/reads.realigned"])

        # sort by positions
        subprocess.call(["samtools",
                         "sort",
                         outDir + "/reads.realigned",
                         "-o",
                         outDir + "/reads.realigned.bam"])

        # index for viewing in IGV
        subprocess.call(["samtools",
                         "index",
                         outDir + "/reads.realigned.bam"])

        # now align to each subtype sequence individually
        for sub in subRefs:
            res = subprocess.check_output(["bowtie2",
                                           "--very-sensitive",
                                           "--no-unal",
                                           "-p",
                                           "4",
                                           "-x",
                                           args[1].rstrip("/") + "/singles/" + sub,
                                           "-U",
                                           outDir + "/reads.fq",
                                           "-S",
                                           outDir + "/" + sub + ".sam"], stderr=subprocess.STDOUT)

            rate = res.split("\n")[-2].split(" ")[0].rstrip("%")
            mult = res.split("\n")[-3].strip().split(" ")[0]
            once = res.split("\n")[-4].strip().split(" ")[0]
            zero = res.split("\n")[-5].strip().split(" ")[0]
            outFP.write(sub + "," + rate + "," + mult + "," + once + "," + zero + "\n")

        outFP.close()


if __name__ == "__main__":
    main(sys.argv[1:])
