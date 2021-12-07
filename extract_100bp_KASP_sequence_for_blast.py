from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import random
import pandas as pd
import os
import subprocess
from collections import defaultdict

with open(sys.argv[1], 'r') as fasta:
	seqs = {str(record.id): str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}


with open(sys.argv[2], 'r') as vcf:
	for line in vcf:
		temp = line.split( )
		chromosome = temp[0]
		start = int(temp[1]) - 51
		end = int(temp[1]) + 50
		print(">" + chromosome + "-" + temp[1])
		print(seqs.get(chromosome)[start:end])