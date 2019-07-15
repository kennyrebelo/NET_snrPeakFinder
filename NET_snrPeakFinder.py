#!/usr/bin/env python


# Import libraries
import argparse
import numpy as np
from os import system



############################################################
#                                                          #
#                           MAIN                           #
#                                                          #
############################################################

# Input arguments
parser = argparse.ArgumentParser(description="""Returns a .bed file per sample
	listing all the single nucleotide peaks discovered.
	run example:
	python NET_snrPeakFinder.py -i NET_Dros_emb_2-3h_S5P_rep1_noSI.bam -g dm6.chrom.sizes""")


parser.add_argument('-i', '--input', dest='input', nargs='+',
                   help="Input One or more .bam file.")

parser.add_argument('-g', '--genome', dest='genome', nargs='?',
                   help="Input genome.chrom.sizes file.")

parser.add_argument('-s', '--strandness', dest='strandness', nargs='?',
                   help="single or paired depending if the data is single or paired-end")



args = parser.parse_args()
bams = args.input
chrom_sizes = args.genome
seq_strategy = args.strandness


for bam in bams:

	# Defining output .bed filename
	filename="Peaks_"+str(bam).split("/")[-1].rstrip(".bam")+".bed"
	output=open(filename,"w")



	# Create .bam files only containing only forward or reverse reads
	if seq_strategy == "single":
		system("samtools view -@ 40 -F 16 -b "+bam+" > F_"+filename+".bam")
		system("samtools view -@ 40 -f 16 -b "+bam+" > R_"+filename+".bam")
		system("samtools sort -@ 40 -o F_"+filename+"_sorted.bam F_"+filename+".bam")
		system("samtools sort -@ 40 -o R_"+filename+"_sorted.bam R_"+filename+".bam")
		system("samtools index F_"+filename+"_sorted.bam")
		system("samtools index R_"+filename+"_sorted.bam")
		system("rm -f F_"+filename+".bam")
		system("rm -f R_"+filename+".bam")
	elif seq_strategy == "paired":
		system("samtools view -@ 40 -f 0x63 -b "+bam+" > F_"+filename+".bam")
		system("samtools view -@ 40 -f 0x53 -b "+bam+" > R_"+filename+".bam")
		system("samtools sort -@ 40 -o F_"+filename+"_sorted.bam F_"+filename+".bam")
		system("samtools sort -@ 40 -o R_"+filename+"_sorted.bam R_"+filename+".bam")
		system("samtools index F_"+filename+"_sorted.bam")
		system("samtools index R_"+filename+"_sorted.bam")
		system("rm -f F_"+filename+".bam")
		system("rm -f R_"+filename+".bam")
	else:
		print "error : select single or paired for the -s argument"


	# Check read distribution over the genome
	system("samtools depth -b "+chrom_sizes+" F_"+filename+"_sorted.bam | awk '$3 > 3' > ntsOver4reads_F.txt")
	system("samtools depth -b "+chrom_sizes+" R_"+filename+"_sorted.bam | awk '$3 > 3' > ntsOver4reads_R.txt")


	nts_F=[]
	with open("ntsOver4reads_F.txt") as file:
		for line in file:
			nts_F.append(line.strip("\n\r").split("\t"))

	nts_R=[]
	with open("ntsOver4reads_R.txt") as file:
		for line in file:
			nts_R.append(line.strip("\n\r").split("\t"))

	nts=map(None, nts_F, nts_R)

	# 
	for nt in nts:


		if nt[0] is not None:

			# Forward reads
			peakPlus=nt[0] #even though they are variables with the same information, it will be useful later
			peakplus_reads=float(nt[0][-1])

			flanku_F=nt[0][0]+":"+str(int(nt[0][1])-100)+"-"+str(int(nt[0][1])-1) #upstream region, considering + strandness (downstream region if - strandness)
			flankd_F=nt[0][0]+":"+str(int(nt[0][1])+1)+"-"+str(int(nt[0][1])+100) #downstream region, considering + strandness (upstream region if - strandness)

			system("samtools depth -aa -r "+flanku_F+" F_"+filename+"_sorted.bam > flanku_F_temp_cov.txt")
			system("samtools depth -aa -r "+flankd_F+" F_"+filename+"_sorted.bam > flankd_F_temp_cov.txt")
			system("cat flanku_F_temp_cov.txt flankd_F_temp_cov.txt > flank_F_temp_cov.txt")

			flankplus_reads=[]

			with open("flank_F_temp_cov.txt") as file:
				for line in file:
					flankplus_reads.append(float(line.strip("\n\r").split("\t")[-1]))

			if peakplus_reads > (3*np.std(flankplus_reads))+np.mean(flankplus_reads):
				peak_coords=peakPlus[0]+"\t"+str(int(peakPlus[1])-1)+"\t"+peakPlus[1]+"\t"+"."+"\t"+"."+"\t"+"+"
				output.write(peak_coords+"\n")


			

		if nt[1] is not None:

			# Reverse reads
			peakMinus=nt[1]
			peakminus_reads=float(nt[1][-1])

			flanku_R=nt[1][0]+":"+str(int(nt[1][1])-100)+"-"+str(int(nt[1][1])-1) #upstream region, considering + strandness (downstream region if - strandness)
			flankd_R=nt[1][0]+":"+str(int(nt[1][1])+1)+"-"+str(int(nt[1][1])+100) #downstream region, considering + strandness (upstream region if - strandness)

			system("samtools depth -aa -r "+flanku_R+" R_"+filename+"_sorted.bam > flanku_R_temp_cov.txt")
			system("samtools depth -aa -r "+flankd_R+" R_"+filename+"_sorted.bam > flankd_R_temp_cov.txt")
			system("cat flanku_R_temp_cov.txt flankd_R_temp_cov.txt > flank_R_temp_cov.txt")

			flankminus_reads=[]

			with open("flank_R_temp_cov.txt") as file:
				for line in file:
					flankminus_reads.append(float(line.strip("\n\r").split("\t")[-1]))

			if peakminus_reads > (3*np.std(flankminus_reads))+np.mean(flankminus_reads):
				peak_coords=peakMinus[0]+"\t"+str(int(peakMinus[1])-1)+"\t"+peakMinus[1]+"\t"+"."+"\t"+"."+"\t"+"-"
				output.write(peak_coords+"\n")






	
	output.close()

system("rm -f nts*")
system("rm -f flank*")
system("rm -f F_*")
system("rm -f R_*")

print "Finished!"
