import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--megaVCF', required=True, type=str,help='path to giant VCF to be chunked')
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-c', '--chunkSize', type=int, help='number of samples included in file', default=1000)
parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)
parser.add_argument('-sd', '--scriptsDirectory', type=str, required=True, help='path to directory where scripts are')


#parser.add_argument('-o', '--outFile', required=True)

'''
THINGS TO NOTICE:
Unix `cut` command 1-indexes file columns

I will be using linux's column ordering system to track columns that are processed 

if one wants to see if a sample ('column') has been processed, it will be the same in the processed file name and the orignial linux ordering

file names are inclusive of all numbers 
'''

args = parser.parse_args()
vcf = args.megaVCF
chunk = args.chunkSize
wd = args.workingDirectory

#makes sure input path wont cause error
if wd[-1] != '/':
    print('add a final slash')
    wd = wd+'/'
sd = args.scriptsDirectory
if sd[-1] != '/':
    print('add a final slash')
    sd = sd+'/'

t = args.threads

#open file to indentify number of samples 
with gzip.open(vcf, 'rt') as v:
    for line in v:
        #if the line is not part of the heading
        if not line.startswith('##'):

            #if the line contains the column names 
            #note that this assumes 9 meta data columns, if you need to account for more/less change this 
            if line.startswith('#'):
                line = line.strip().split('\t')
                #numSamps = len(line[9:])
                lenRow = len(line)

                #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                print('lenRow',lenRow)
                break
        
        #else:
        #    print(line)


#old code, definitely too complex
'''
if numSamps % chunk != 0:
    remainder = (numSamps % chunk)
    print('remainder', remainder)
    #dk if i need upper
    upper = numSamps + remainder
else:
    upper = numSamps
#columns are 1-indexed 

print('ipper', upper)
print('num samps',numSamps, '\n')
'''


#note this assumes 9 meta cols, will need to change here as well 
#this is all 1-indexing dont get confused
for i in range(10,lenRow,chunk):
    #print('i',i)


    #split up the total number of samps (all columns - metadata cols) in as equal chunks as you can
    #this means that the starting and ending chunks will have fewer samples than the middle chunks

    #for the last set of samples
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
    
    print('wd', wd)
    firstOut = f'{start}-{end}samps.vcf'
    secondOut = f'{start}-{end}samps.filtered.vcf'
    thirdOut = f'{start}-{end}samps.filtered.counted.vcf'
    fourthOut = f'{start}-{end}samps.filtered.counted.processed.vcf'
    

    #print('first out', firstOut)
    #print('second out', secondOut)
    #print('\n')
    
    #create the chunk file 
    #immediately bgzip it
    #filter w bcftools
    #gzip
    #process_VCF creates allele counts and deletes useless info
    #gzip
    #run process_megaVCF on it, deletes useless info and processes indels 

    #start never = 1... delete?
    if start == 1:
        os.system(f'gzip -dc {vcf} | cut -f1-{end} > {wd}{firstOut}')
        #test without indexing, then index if it goes badly
        os.system(f'bgzip -@ {t} -f {wd}{firstOut}')
        os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}{secondOut} {wd}{firstOut}.gz")

        #def delete first out after filtering!!!!
        #os.system(f'rm {wd}{firstOut}.gz')
        
        os.system(f'gzip -f {wd}{secondOut}')
        os.system(f"python3 {sd}process_VCF.py -v {wd}{secondOut}.gz -wd {wd} -o {thirdOut}")
        #os.system(f'rm {wd}{secondOut}.gz')
        os.system(f'gzip -f {wd}{thirdOut}')

        os.system(f"python3 {sd}current-process-vcf.py -v {wd}{thirdOut}.gz -o {wd}{fourthOut}")

        #os.system(f"rm " first out)
        #os.system(f"gzip secondout")
    else:
        os.system(f'gzip -dc {vcf} | cut -f1-9,{start}-{end} > {wd}{firstOut}')
        os.system(f'bgzip -@ {t} -f {wd}{firstOut}')

        #make sure this file stays here or figure out a better way to do this 
        os.system(f"bcftools annotate -x '^FORMAT/GT' -h /more_storage/lily/tb_tree/parse_vcf/add-header -O v -o {wd}{secondOut} {wd}{firstOut}.gz")

        #os.system(f'rm {wd}{firstOut}.gz')
        os.system(f'gzip -f {wd}{secondOut}')
        os.system(f"python3 {sd}process_VCF.py -v {wd}{secondOut}.gz -wd {wd} -o {thirdOut}")
        #os.system(f'rm {wd}{secondOut}.gz')
        os.system(f'gzip -f {wd}{thirdOut}')


        os.system(f"python3 {sd}current-process-vcf.py -v {wd}{thirdOut}.gz -o {wd}{fourthOut}")
        #os.system(f"rm " first out)
        #os.system(f"gzip secondout")
        
    












#old code that is too complex


"""
for i in range(1,upper,chunk):
    print('i',i)
    count += 1
    #this means that the starting and ending chunks will have fewer samples than the middle chunks
    # maybe change this later if it matters 
    #start = i   
    if i == 1:
        start = 10
        end = 10 + chunk 

    elif i+chunk > numSamps:
        print('subset')
        start = i+10
        end = upper+1
        
    else:
        print('most of these?')
        start = i+10
        end = i+chunk-1+10
        
    print(start, end)
    #create the chunk file 
    #immediately bgzip it
    #filter w bcftools
    #switch to gzip?
    #run process_megaVCF on it 
    print('wd', wd)
    firstOut = f'{start}-{end}samps.vcf'
    secondOut = f'{start}-{end}samps.filtered.vcf'
    thirdOut = f'{start}-{end}samps.filtered.counted.vcf'
    #firstOut = f'{start}-{end}samps.vcf'
    #secondOut = f'{start}-{end}samps.filtered.vcf.gz'
    print('first out', firstOut)
    print('second out', secondOut)
    
    if start == 1:
        os.system(f'gzip -dc {vcf} | cut -f1-{end} > {wd}{firstOut}')
        #test without indexing, then index if it goes badly
        os.system(f'bgzip -@ {t} {wd}{firstOut}')
        os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}{secondOut} {wd}{firstOut}.gz")
        #os.system(f'rm {wd}{firstOut}.gz')
        os.system(f'gzip {wd}{secondOut}')
        os.system(f"python3 {sd}process_VCF.py -v {secondOut}.gz -wd {wd} -o {thirdOut}")
        #os.system(f'rm {wd}{secondOut}.gz')
        os.system(f'gzip {wd}{thirdOut}')

        #os.system(f"rm " first out)
        #os.system(f"gzip secondout")
    else:
        os.system(f'gzip -dc {vcf} | cut -f1-9,{start}-{end} > {wd}{firstOut}')
        os.system(f'bgzip -@ {t} {wd}{firstOut}')

        #make sure this file stays here or figure out a better way to do this 
        os.system(f"bcftools annotate -x '^FORMAT/GT' -h /more_storage/lily/tb_tree/parse_vcf/add-header -O v -o {wd}{secondOut} {wd}{firstOut}.gz")

        #os.system(f'rm {wd}{firstOut}.gz')
        os.system(f'bgzip {wd}{secondOut}')
        os.system(f"python3 {sd}process_VCF.py -v {secondOut}.gz -wd {wd} -o {thirdOut}")
        #os.system(f'rm {wd}{secondOut}.gz')
        os.system(f'gzip {wd}{thirdOut}')
        #os.system(f"rm " first out)
        #os.system(f"gzip secondout")
        """
    
    
    

    #failproof
    #print(start,end)
    #if count == 1:
    #    print('coubt', count)
    #    break
    #print('no break')



      





#os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {out} {vcf}")
