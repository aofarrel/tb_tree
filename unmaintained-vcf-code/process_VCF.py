import sys
import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-wd', '--working_directory', required=True, help='path to working directory')
parser.add_argument('-v', '--vcf_file', required=True, help='VCF input')
parser.add_argument('-o', '--out_file', required=True, help='processed VCF output')
args = parser.parse_args()

vcf=args.vcf_file
print(vcf)
wd = args.working_directory
if wd[-1] != '/':
    #print('add a final slash')
    wd = wd+'/'
out=args.out_file

'''
THIS PROGRAM IS NO LONGER USED IN THE PROCESSING OF TB VCFS FOR USHER PLACEMENT. 
KEPT FOR RECORDS BUT NO LONGER USED 
'''
#this program will take a chunked vcfg post-filtering and preprocessing
#note this code also assumes 9 cols of meta data 
with gzip.open(vcf, 'rt') as v:
    with open(out, 'w') as o:
        kept = 0
        removed = 0
        for line in v:
            #skip through header
            if not line.startswith('##'):
                #find row with sample names
                if line.startswith('#'):
                    #this should write to file in correct order
                    o.write(line)
                    line = line.strip().split('\t')
                    #keep track of sample names (note that these are 0 indexed)
                    sampNames = line

                    #keep an ongoing count of missing data for each sample
                    #this is 0-indexed!!!! needs to be 1-indexed for unix `cut`
                    samples = {i:0 for i in range(9,len(line))}

                #all position data goes here
                #process each line for allele counts
                else:
                    line = line.strip().split('\t')
                    
                    #delete this after testing 
                    #if int(line[1]) > 3000:
                    #    print('done')
                    #    break
                    
                    #figure out how many alt alleles are at this position
                    alleles = line[4].split(',')

                    #track the allele counts for the position
                    counts = {(i+1):0 for i in range(len(alleles))}

                    #go through all samples for the position to determine counts
                    #use include to track if a line should be deleted from the file 
                    include = False
                    an = 0
                    for i in range(9, len(line)):
                        al = line[i].split('/')
                        
                        #there shouldnt be any heterozygous alleles
                        assert al[0] == al[1]
                        
                        #if there is missing data, need to include the line in the VCF
                        #if there is missing data note which sample it's in 
                        if al[0] == '.':
                            #print('missing', line[1])
                            include = True
                            samples[i] += 1
                        
                        #keep track of all alt alleles that have a count > 0
                        else:
                            an += 1
                            if int(al[0]) > 0:
                                counts[int(al[0])] += 1

                    #keys must be ints to sort correctly
                    keys = sorted(counts.keys())
                    ac = ''

                    for c in keys:
                        #conditional for formatting 
                        if len(ac) == 0:
                            ac += str(counts[c])
                        else:
                            ac += ','+ str(counts[c])

                        #if any sample has a variant at this position, the line must be kept in the file
                        if counts[c] > 0:
                            include = True
                            #print('variant', line[1])
                        
                    #print('ac', ac)
                    #if a variant or missing data was found, include = True and the line is kept
                    if include == True:
                        kept += 1
                        #print(ac,an)
                        #update AC and AN and write line to file 
                        line[7] = 'AC='+ ac + ';' + 'AN=' + str(an)
                        o.write('\t'.join(line)+'\n')  
                    else:
                        removed += 1
                        #print('exclude', line[1])
            #write header to new file
            else:
                o.write(line)


#iterate through samples to determine which (if any)should be removed 
#print('kept', kept) 
#print('removed', removed ) 
total = kept+removed
#print('total', total ) 
cols = sorted(list(samples.keys()))
#print('samp names', sampNames)
#print('cols',cols)
deleteCols = []
#to check for shitty samples
for c in cols:
    #print(out)
    #print('c', c, samples[c])
    #change to 1-index
    if samples[c]/total > .05:
        #c is the actual column to delete
        c = c+1 
        deleteCols.append(c)

#remove samples from file
#record all samples removed in file (can combine these files later)
#print('to delete',deleteCols)
#print('end cols', cols[-1])
if len(deleteCols) > 0:
    #print('MODIFYING FILE')
    with open(f'{wd}deletedsamplesin{out}.txt','w') as ds:
        #if there is only one sample that needs to be deleted, cut is a faster method
        if len(deleteCols) == 1:
            ds.write(f'{sampNames[deleteCols[0]-1]}\n')
            #assumes 9 columns of metadata and 1 indexing for cols
            #file chunks always have first sample at column 10 (1 ind)
            if deleteCols[0] == 10:
                #print(f'cut -f1-9,11-{cols[-1]+1} {wd}{out} > {wd}current.mod')
                os.system(f'cut -f1-9,11-{cols[-1]+1} {wd}{out} > {wd}current.mod')
                #print('overwrite')
                os.system(f'mv current.mod {out}') 
            
            #if the column to be deleted is equivlant to last column in cols (note cols is still 0-indexed)
            elif deleteCols[0] == cols[-1]+1:
                #print(f'cut -f1-{deleteCols[0]-1} {wd}{out} > {wd}current.mod')
                os.system(f'cut -f1-{deleteCols[0]-1} {wd}{out} > {wd}current.mod')
                #print('overwrite')
                os.system(f'mv {wd}current.mod {wd}{out}')
            #if the col to be deleted is somewhere in the middle
            else:
                #print(f'cut -f1-{deleteCols[0]-1},{deleteCols[0]+1}-{cols[-1]+1} {wd}{out} > {wd}current.mod')
                os.system(f'cut -f1-{deleteCols[0]-1},{deleteCols[0]+1}-{cols[-1]+1} {wd}{out} > {wd}current.mod')
                #print('overwrite')
                os.system(f'mv {wd}current.mod {wd}{out}')

        #if there is more than one column to be deleted, cut function will cause too many probs
        #record deleted samps in file and use bcftools
        else:
            #print('make delete file')
            for d in deleteCols:
                #print(d)
                ds.write(f'{sampNames[d-1]}\n')
    
    #for all files with more than 1 sample deletion
    if len(deleteCols) > 1:
        #print(f'{wd}{out}')
        #print(f'bcftools view -S ^{wd}deletedsamplesin{out}.txt -O v -o {wd}current.mod -I {wd}{out}')
        os.system(f'bcftools view -S ^{wd}deletedsamplesin{out}.txt -O v -o {wd}current.mod -I {wd}{out}')
        #os.system(f'bcftools view -S ^{wd}deletedsamplesin{out}.txt -O v -o {wd}current.mod {wd}{out}')
        #print('overwrite')
        os.system(f'cp {wd}current.mod {wd}{out}currentmod')
        os.system(f'mv {wd}current.mod {wd}{out}')
        #print('overwrite')
        #os.system(f'mv {wd}current.mod {wd}{out}')
        #break
            
            
print('done w process_VCF', 'total lines:', total, 'lines kept:', kept, 'lines removed:', removed)        

                
            
            
            
                
                
            
            
                   
                    

                    