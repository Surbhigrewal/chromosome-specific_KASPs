from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

with open(sys.argv[1], 'r') as fasta:
	seqs = {str(record.id): str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}

with open(sys.argv[2], 'r') as vcf:
	ref_alt_allele_dict = {}
	for line in vcf:
		temp = line.split( )
		ref_alt_allele_dict[temp[0]+"-"+temp[1]] = (temp[3], temp[4])

print("\t".join(["SNP_ID", "Chr", "Pos", "REF", "ALT", "SNP_Seq"]))
with open(sys.argv[3], 'r') as id_file:
	for line in id_file:
		chromosome = line.split("-")[0]
		snp_site = line.rstrip( ).split("-")[1]
		alleles = ref_alt_allele_dict.get(line.rstrip( ))
		ref_allele = alleles[0]
		alt_allele = alleles[1]
		seq = seqs.get(line.rstrip( ))
		print(line.rstrip( ) + "\t" + chromosome.split("chr")[1] + "\t" + snp_site + "\t" + ref_allele + "\t" + alt_allele + "\t" + seq[0:50] + "[" + ref_allele + "/" + alt_allele + "]" + seq[51:])
