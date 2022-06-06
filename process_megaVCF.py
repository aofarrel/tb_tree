import sys
import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf_file', required=True, help='VCF input')
parser.add_argument('-o', '--out_file', required=True, help='processed VCF output')
args = parser.parse_args()
vcf=args.vcf_file
out=args.out_file

#run once to extract the right filters 
#os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {out} {vcf}")
#
#os.system(f"zcat {vcf} | gzip -c > {vcf}.nob.gz")
#the file in this should be the file created from the system command above
with open(vcf) as v:
    for line in v:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            print(line)
            alleles = line[4].split(',')
            counts = {i+1 for i in range(len(alleles))}
            print('counts', counts)
            for i in range(9, len(line)):
                al = line[i].split('/')
                if int(al[0]) > 0:
                    print(line[i])
                    