import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--samples_file', required=True)
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-o', '--srr_file', required=True)

args = parser.parse_args()
samps=args.samples_file
srrs=args.srr_file

with open(samps) as f:
    count = 0 
    for line in f:
        line = line.strip().split('.')
        samp = line[5]
        os.system(f'ffq {samp} -o {samp}.json')
        os.system(f"cat {samp}.json | grep accession | perl -pe 's/\"//g' | perl -pe 's/,//g' | perl -pe 's/accession/\naccession/g' | grep 'SRR' | cut -d ' ' -f2 | uniq > {samp}accession")
        
    os.system(f'cat *accession > SRAaccessions')


    


    
        


