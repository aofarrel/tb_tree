import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--megaVCF', required=True, type=str,help='path to giant VCF to be chunked')
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-c', '--chunkSize', type=int, help='number of samples included in file', default=1000)
parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)


#parser.add_argument('-o', '--outFile', required=True)
args = parser.parse_args()
vcf = args.megaVCF
chunk = args.chunkSize
wd = args.workingDirectory
t = args.threads
#outFile=args.outFile

with gzip.open(vcf, 'rt') as v:
    for line in v:
        #print(line)
        if not line.startswith('##'):

            if line.startswith('#'):
                line = line.strip().split('\t')
                print(line[:10])
                numSamps = len(line[9:])
                break
        else:
            print(line)

print(numSamps)
count = 0
for i in range(0,numSamps,chunk):
    count += 1
    #this means that the starting and ending chunks will have fewer samples than the middle chunks
    # maybe change this later if it matters 
    start = i    
    if i+chunk > numSamps:
        end = numSamps
    else:
        end = i+chunk
    
    #create the chunk file 
    #immediately bgzip it
    #filter w bcftools
    #switch to gzip?
    #run process_megaVCF on it 
    firstOut = f'{start}-{end}samps.vcf'
    secondOut = f'{start}-{end}samps.filtered.vcf.gz'
    
    if start == 0:
        os.system(f'gzip -dc {vcf} | cut -f1-{end} > {wd}/{firstOut}')
        #test without indexing, then index if it goes badly
        os.system(f'bgzip -@ {t} {wd}/{firstOut}')
        os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}/{secondOut} {wd}/{firstOut}.gz")
    else:
        os.system(f'gzip -dc {vcf} | cut -f1-9,{start}-{end} > {wd}/{firstOut}')
        os.system(f'bgzip -@ {t} {wd}/{firstOut}')
        os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}/{secondOut} {wd}/{firstOut}.gz")
    

    #failproof
    #print(start,end)
    
    if count == 2:
        print('coubt', count)
        break



        





#os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {out} {vcf}")
