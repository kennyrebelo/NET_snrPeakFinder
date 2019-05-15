# NET_snrPeakFinder

Python script that from specific regions of interest and a single nucleotide resolution NET-seq alignment file returns only the enriched regions, or peaks. Peak defined as in *insert reference*

Described as part of an analysis pipeline in *insert reference*

Parameters:

-i, --input	input .bam files

-g, --genome	genome regions of interest in .bed format

-s, --strategy	single/paired 


Usage example: 

python NET_snrPeakFinder.py -i sample1.bam sample2.bam -g regions.bed -s paired
