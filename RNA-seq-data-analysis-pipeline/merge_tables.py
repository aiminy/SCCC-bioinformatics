#!/usr/bin/python
#This code is adapted from Dave Wheeler
#To run for our data:
#make a file called Sample_10_raw_count.txt that include the following 
#note count file name and sample name are separated by space instead of TAB
#201348189-01_raw_count_0.txt cbpKO_47_sham
#201348190-01_raw_count_0.txt cbpKO_36_TAC
#201348191-01_raw_count_0.txt cbpKO_32_TAC
#201348192-01_raw_count_0.txt WT_12_sham
#201348193-01_raw_count_0.txt WT_21_sham
#201348195-01_raw_count_0.txt WT_4_TAC
#201348196-01_raw_count_0.txt p300KO_67_sham
#201348199-01_raw_count_0.txt p300KO_79_TAC
#201348200-01_raw_count_0.txt cbpKO_33_sham
# then run

#python python ~/Code/RNA-seq-data-analysis-pipeline/merge_tables.py Sample_10_raw_count.txt Output_count_samples10.txt

# you will get file name as "Output_count_samples10.txt"
# this file is used as input for Run_DESeq2.R

#####################################################
#   example.py - a program to ....                  #
#                                                   #
# Author: Dave Wheeler                              #
#                                                   #
# Purpose: merge count tables                       #
#                                                   #
# Usage: python merge_tables guide_file             #
#####################################################
#looks for text guide_file
#that contains files to be merged with column headers (space separted) 
#ie
#file1.counts untreated1
#file2.counts untreated2
#file3.counts treated1
#file4.counts treated2
#this will generated a tab separated table like this
#
#gene untreated1 untreated2 treated1 treated2
#gene1 0 0 0 0 
#gene2 1 0 11 10
#.......
##############################################
import sys
try:
	infile = open(sys.argv[1])
except IndexError:
	print "No guide file provided"
	sys.exit()

outf = sys.argv[2]
	
#make dict of genes with list of counts
#list is ordered so treatments will be preserved.
#genes = {'gene1':[1,2,3,4]}
#header keeps track of treatment order, will be as read from config
col_header = []
genes = {}	

outfile = open(outf,'w')

for line in infile:
	filename,header = line.strip().split(' ')
	try:
		data_f = open(filename)
	except IOError:
		print "%s can't be found?"%filename
		sys.exit()
		
	col_header.append(header)	
	
	#read file and add gene and counts to the dict
	for line in data_f:
		gene,count = line.strip().split('\t')
		if gene not in genes:
			genes[gene] = [count]
		else:
			genes[gene].append(count)
	#important to close file
	data_f.close()	
	
infile.close()

outfile.write('gene\t'+'\t'.join(col_header)+'\n')

for gene in genes:
	data = genes[gene]
	#make sure each treatment has a count for this gene
	#this should catch most errors
	try:
		assert len(data) == len(col_header)
	except AssertionError:
		print "one of the treatment or genes is missing or extra"
		print "data, found the problem here:"
		print gene,data
		print "while %s columns of treatments given" %len(col_header)
		sys.exit()
		
	out_data = gene+'\t'+'\t'.join(data)+'\n'
	outfile.write(out_data)
	
outfile.close()
#print "Merged table is 'merged_counts.txt'"
print 'Merged table is:', outf
