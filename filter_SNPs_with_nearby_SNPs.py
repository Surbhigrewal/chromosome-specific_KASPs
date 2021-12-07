import argparse
parser = argparse.ArgumentParser(description="identify snps in muticum vcf that don't have anthoerh snp closeby")
parser.add_argument("input_vcf", action="store", help="vcf of muticum INCLUDING HET AND NON-PASSES")
parser.add_argument("dist", action="store", help="distance between snps", default= 50)
parser.add_argument("output_vcf", action="store", help="output file", default= 50)

args = parser.parse_args()

f1 = open(args.output_vcf, "w")
listy = []
yn = 0
with open(args.input_vcf) as file:
    first = file.readline()
    pos = int(first.split()[1])
    for line in file:
        plus50 = pos + int(args.dist)
        neg50 = pos - int(args.dist)
        if yn == 1:
            if int(line.split()[1]) < (int(listy[-1].split()[1]) +int(args.dist)):
                del listy[-1]
        if int(line.split()[1]) >= plus50:
            listy.append(line)
            yn = 1
        else: yn = 0
        pos = int(line.split()[1])

for item in listy:
    f1.write(item)