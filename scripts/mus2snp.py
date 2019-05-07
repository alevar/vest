#!/usr/bin/env python

import sys

import pandas as pd

import argparse


def mus2snp(args):
    genomes = dict()

    curGenome = ""

    with open(args.input, "r") as musFP:
        for line in musFP.readlines():
            if line[0] == ">":
                curGenome = line[1:].split()[0]
                genomes.setdefault(curGenome, "")
            else:
                assert not curGenome == "", "empty genome name"
                genomes[curGenome] += line.strip().upper()

    # assert that all sequences are of the same length
    count = 0
    curLen = 0
    for gname in genomes:
        if count == 0:
            curLen = len(genomes[gname])
        else:
            assert len(genomes[gname]) == curLen, "genomes have different lengths"

        count += 1

    # now to actually extract the SNPs and indels

    ref_seq = genomes[args.refname]

    snpList = list()

    id_count = 0

    for gname in genomes:
        cur_seq = genomes[gname]
        if not gname == args.refname:  # make sure we are not performing a self-to-self comparison

            cur_ins = ""  # keeps track of the inserted sequence

            cur_del = 0  # keeps track of the length of the current deletion

            ref_pos = 0  # reset reference position

            for pos in range(curLen):  # iterate over all positions
                if not ref_seq[pos] == "-":  # increment the reference position for each base that is present
                    ref_pos += 1

                id_count += 1
                if not ref_seq[pos] == "-" and not len(cur_ins) == 0:  # all insertions have been processed
                    cur_ins = cur_ins.replace("-", "")
                    if not len(cur_ins) == 0:  # does an actual insertion exist
                        if args.ins:
                            snpList.append([args.refname + "__" + gname + "__" + str(id_count),
                                            "insertion",
                                            args.refname,
                                            ref_pos,
                                            cur_ins])
                        cur_ins = ""

                if not cur_seq[pos] == "-" and not cur_del == 0:  # all deletions have been processed
                    if args.dele:
                        snpList.append([args.refname + "__" + gname + "__" + str(id_count),
                                        "deletion",
                                        args.refname,
                                        ref_pos - cur_del,
                                        cur_del])
                    cur_del = 0

                if ref_seq[pos] == "-":  # there is an insertion
                    cur_ins += cur_seq[pos]

                elif cur_seq[pos] == "-":  # there is a deletion
                    cur_del += 1

                elif not ref_seq[pos] == cur_seq[pos]:  # there is a SNP
                    if not args.nosnp:
                        snpList.append([args.refname + "__" + gname + "__" + str(id_count),
                                        "single",
                                        args.refname,
                                        ref_pos,
                                        cur_seq[pos]])
                else:  # there is a match
                    continue

    snpDF = pd.DataFrame(snpList, columns=["name", "type", "chr", "pos", "last"])
    snpDF.drop_duplicates(["chr", "pos", "last"], keep='first', inplace=True)
    snpDF.sort_values(by="pos", inplace=True)
    snpDF.reset_index(drop=True, inplace=True)
    snpDF.to_csv(args.output, sep="\t", header=False, index=False)


def main(argv):
    parser = argparse.ArgumentParser(description='''Help Page''')

    # ===========================================
    # ===================BUILD===================
    # ===========================================
    parser.add_argument('--input',
                        required=True,
                        type=str,
                        help="multiple sequence alignment in mus format")
    parser.add_argument('--refname',
                        required=True,
                        type=str,
                        help="path to the genome sequence corresponding to the annotation")
    parser.add_argument('--output',
                        required=True,
                        type=str,
                        help="output file path")
    parser.add_argument('--nosnp',
                        required=False,
                        action='store_true',
                        help="do not include SNPs")
    parser.add_argument('--ins',
                        required=False,
                        action='store_true',
                        help="include insertions")
    parser.add_argument('--dele',
                        required=False,
                        action='store_true',
                        help="include deletions")

    parser.set_defaults(func=mus2snp)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
