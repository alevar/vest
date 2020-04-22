# VEST
VEST is a fast and accurate method for increasing mapping rate of RNA and DNA short reads sequenced 
from hypermutated genomes by performing graph-assisted realignment of mappings produced by conventional 
aligners (Bowtie2,HISAT2,STAR,BWA,MiniMap2,etc). The method is primarily developed to assist in experiments 
where sequencing depth is low and assembly of the genome may not be feasible (single-cell experiments, 
metagenomic analysis, etc).

## Installation
In most cases Vest should be isntalled from source with the following steps:
1. cmake -DCMAKE_BUILD_TYPE=Release .
2. make
3. make install

## VEST Help Page

```asm
Vest
Modes:
	build - build the graph from the MSA
	realign - map alignments onto the MSA and infer consensus
	inspect - inspect the graph
	build-refind - build additional indices for the refind subcommand used to identify recombination events
	refind - build recombination map for a given genome
```

## Preparing data for indexing with VEST

Vest requires a Multiple Sequence Alignment (MSA) as input to build the population graph. MSA can be built 
using a variety of methods such as MUSCLE and MAFFT which output a FASTA-formatted alignment.

#### Building Index

```asm
vest_build Help Page

vest_build
-i -o 
Arguments:
	i/--mus	path to the multiple sequence alignment
	o/--out	path and basename for the output database
```

###### Preparing Alignments
Vest does not perform alignment of the reads, allowing flexibility of chosing whichever option works best 
for a given experiment. As such, one must build a separate index to align reads against as per documentation 
of the chose method. The index must be build using the same set of sequences (or a subset) used to build 
the index for Vest. One may obtain the fasta file for the reference population index by using `vest inspect`.

Reads can then be aligned using the method of choosing to prepare the input for Vest. it is required, however, 
that the input alignment is SAM/BAM/CRAM formatted, as other formats are currently not supported.
 
#### Realignment of Reads

Once the index and the input alignment is obtained, one may perform the realignment with vest by using 
`vest realign` command.

```asm
vest_realign Help Page

vest_realign
-o -x [-a -b -i -k -n -p -t ]
Arguments:
	a/--ann	annotation of genomic features with respect to one of the genomes
	b/--bed	bed-formatted interval to be converted to the reference sequence
	i/--input	path to the input alignment
	k/--keep	Keep tmp data
	n/--gname	Name of one of the sequences in the pan genome to be used when resolving a gap in the assembly. If not provided, the default behavior is to use the sequence supported by the reads on both sides of the junction
	o/--output	base name for the output files
	p/--pname	Name of the sequence in the pan genome to which the alignments will be projected
	t/--tmp_dir	directory for storing temporary data
	x/--db	database build with vest build from the ultiple sequence alignment
```

Several options are available here to control how alignments are fitted onto the graph, and 
what additional information can be inferred in the process.

Often, full coverage will not be possible to create a 100% contiguous assembly. Several strategies are 
currently available for filling the gaps:
1. Using predefined reference genome. This can be used to enforce that gaps are filled with unused nucleotides 
from a provided genome (such as HXB2 for example)
2. By default, Vest will scan 3' and 5' covered regions of the alignment and fill the gap with the most abundant 
genome at those positions

Additionally, if an anotation in GTF/GFF format exists for one of the genomes in the index, it can be 
fitted to the inferred consensus sequence in the process of realignment. The annotation file can be provided 
with `-a` option.