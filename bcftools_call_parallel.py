import argparse
import multiprocessing
import os
parser = argparse.ArgumentParser(description="Process command input")
parser.add_argument("bam", action="store", help="bam with csi index")
parser.add_argument("prefix", action="store", default=10, help="prefix of snp calling run")
args = parser.parse_args()

def worker(com):
    """runs the command item in os.system"""
    os.system(com)


chroms = ["chr1A", "chr1B", "chr1D", "chr2A", "chr2B", "chr2D", "chr3A", "chr3B", "chr3D", "chr4A", "chr4B", "chr4D", "chr5A", "chr5B", "chr5D", "chr6A", "chr6B", "chr6D", "chr7A", "chr7B", "chr7D", "chrUn"]

commands = []
for item in chroms:
    commands.append("samtools mpileup -v --output-tags AD,DP -r {0}  -f iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta {1} | bcftools call -mv -Ov -o {2}_variants_{0}.raw.vcf"
                .format(item, args.bam, args.prefix))

if __name__ == '__main__':
    jobs = []
    for item in commands:
        p = multiprocessing.Process(target=worker, args=(item,))
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()
