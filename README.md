# NET_snrPeakFinder

Python script that requires samtools installed in run environment. Developed and tested with [SAMtools](http://samtools.sourceforge.net/) v1.7.

From a set of genomic regions of interest (.bed) it will return the nucleotide regions (.bed) that were enriched for NET-seq signal (i.e., nucleotides with a NET-seq peak).


Parameters:

-i, --input	input .bam files

-g, --genome	genome regions of interest in .bed format

-s, --strategy	single/paired 

Designed for NET-seq data and described as part of an analysis pipeline in *insert reference*

Usage example: 

python NET_snrPeakFinder.py -i mNET_Long_S5P_rep1_SNR.bam -g gene_list.bed -s paired


### It operates in the following manner:

Starts by separating forward strand aligning reads from reverse strand aligning reads. SAM flags used will be appropriate for single-end or paired-end reads.

Read coverage will be inquired for the input .bed regions (for both strands) and only the single nucleotide regions that have at least 4 reads will be saved into a file. One “F” file for forward strand aligning reads and one “R” file for reverse strand aligning reads. In this specific case since we are using samtools depth on alignment files (.bam) with single nucleotide sized reads the output for the depth command will be the the single nucleotide positions (1-based) of the alignment file for the regions provided as input. From which we only save the positions with 4 reads or more obtaining thus a list of potential candidates for a peak.

Then for each possible peak-candidate we check coverage values for the surrounding 200 nucleotides and calculate the mean and standard deviation values for this flanking region. If the read number in the candidate is bigger than the mean plus 3 times the standard deviation it is considered a peak. This is done for both strands.

```
Peak > mean + (3*standard deviation)
```
