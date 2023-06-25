#!/usr/bin/python
import gzip
import sys
from Bio import SeqIO
import argparse

'''
filterFastqNReads.py
FUNCTION:
This script will take the next generation sequencing reads in fastq format. 
          Then it will filter out the sequence which has unknown base pair Ns. 
          The reads without N will be output to the result file while those with Ns will be discarded. 
          This will retrieve the perfect reads without missing data.
          Statistcal summary will be generated for each input seq file

HELP: python3 Fastq_stat.py -h [--help]

USAGE: 

python3 Fastq_stat.py -i input_fastq_file -o out_statistical_result_file

Example on testing data:

python3 Fastq_stat.py -i testdata_illumina_1.fq -o fastq.stat.txt

Path information could be added before the input and output file name if files are not in current directory.

Version 2021-June-15th

'''
__author__ ="Xuewen Wang"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",required=True, help="The input file containing reads in fastq format")
parser.add_argument("-o","--output", default='fastq.stat.txt', help="Statistical of read file")
args=parser.parse_args()
ind=args.input
otd=args.output


lengths = []
sumlen = 0
sumraw=0
ct=0
ctgood=0
subword="N"

# outf = open(otd, 'w')
with gzip.open(ind,"rt") as IN:
 for record in SeqIO.parse(IN, "fastq"):
     ct +=1
     sumraw += len(record.seq)
     if subword in record.seq:
          next
     else:     
          # SeqIO.write(record, outf, "fastq")
          sumlen += len(record.seq)
          ctgood +=1


#average length
lenmean=float("{:.1f}".format(sumraw/ct))


# Report stats
# outf.write("Sequence file:\t"+ind)
print ("Sequence file:\t", ind)
print("Total number of input sequences/reads:\t", ct)
print("Total length (bp) of input sequences/reads:\t", sumraw)
print("Mean length (bp) of input sequences/reads:\t", lenmean)
print("Total number of sequences/reads after filtering:\t", ctgood)
print("Total length (bp) of filtered sequences/reads:\t", sumlen)
#print("Result:\t",otd)

