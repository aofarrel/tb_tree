import argparse
import gzip as gz
import os

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcfFile', required=True)
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-c', '--chunk', type=int,required=True)


args = parser.parse_args()
vcf=args.vcfFile
chunk=args.chunk

wd = args.workingDirectory


'''
THIS SCRIPT IS VERY SIMPLE AND PURELY DESIGNED TO REDUCE THE SIZE OF THE INPUT VCF TO REDUCE PROCESSING TIME. THIS SCRIPT IS NOT WELL DEVELOPED AND SHOULD NOT BE 
RUN WITHOUT TESTING AND DETERMINING HOW MUCH STORAGE SPACE IS NEEDED. THIS SCRIPT IS NOT OPTIMIZED OR PARALLELIZED AND WILL LIKELY ONLY BE IF ITS WORTH THE TIME IT SAVES. 
IT IS IMPORTANT TO NOTICE THAT RUNNING THIS SCRIPT ON THE SAME FILE WILL CAUSE A SIGNIFICANT SLOWDOWN AS SCRIPTS THAT ACT ON THE SAME FILE WILL RUN INTO IO LIMITATIONS.

DEPENDING ON THE SIZE OF THE VCF, 1000 SAMPLES IS A GOOD CHUNK SIZE AS IT WILL CONSIDERABLY REDUCE THE FILE SIZE WITHOUT RESULTING IN AN EXTREME NUMBER OF FILES 
'''

#makes sure input path wont cause error
if wd[-1] != '/':
    print('add a final slash')
    wd = wd+'/'

samples = {}
with gz.open(vcf, 'rt') as v:
    for line in v:
        #if the line is not part of the heading
        if not line.startswith('##'):

            #if the line contains the column names 
            #note that this assumes 9 meta data columns, if you need to account for more/less change this 
            if line.startswith('#'):
                line = line.strip().split('\t')
                print(line)
                for i in range(9,len(line[9:])):
                    #print(i)
                    samples[i] = line[i]
                #print(samples)
                #numSamps = len(line[9:])
                lenRow = len(line)

                #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                #print('lenRow',lenRow)
                break
        
        #else:
        #    print(line)


#note this assumes 9 meta cols, will need to change here as well 
#this is all 1-indexing dont get confused

#change this loop after testing 

#for i in range(10,lenRow):
for i in range(10,lenRow, chunk):
    print(i)
    print(i+chunk)

    if i+chunk > lenRow:
        print('subset')
        start = i
        end = lenRow
    #for the first set of samples (assumes 9 cols of metadata)
    elif i == 10:
        start = i
        end = i + chunk 
    #all middle cols
    else:
        start = i+1
        end = i + chunk 

    print(start, '-', end)
    
    #need to cut these into chunks of even size of lenrow 
    os.system(f'gzip -dc {vcf} | cut -f1-9,{start}-{end} > {wd}col{start}-{end}.vcf')
    os.system(f'bgzip -f {wd}col{start}-{end}.vcf')
    
    
    


