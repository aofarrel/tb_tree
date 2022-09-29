import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', required=True, type=str,help='path to giant VCF to be chunked')
parser.add_argument('-d', '--work_dir', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)
#parser.add_argument('-sd', '--scripts-dir', type=str, required=True, help='path to directory where scripts are')


'''
THINGS TO NOTICE:
Unix `cut` command 1-indexes file columns

I will be using linux's column ordering system to track columns that are processed 

if one wants to see if a sample ('column') has been processed, it will be the same in the processed file name and the orignial linux ordering

file names are inclusive of all numbers 

diff file lists 1st position and the full lenght of the string w this char, this may need to change becasue the lenght value is inclusive of the first sight 
'''

args = parser.parse_args()
vcf = args.vcf
wd = args.work_dir

#makes sure input path wont cause error
if wd[-1] != '/':
    print('add a final slash')
    wd = wd+'/'
#sd = args.scriptsDirectory
#if sd[-1] != '/':
#    print('add a final slash')
#    sd = sd+'/'

t = args.threads


def find_snps(line):
    ref = line[3]
    alt = line[4]
    end = len(ref)
    snps = []
    for i in range(end):
        if ref[i] != alt[i]:
            snps.append(i)
    
        
    print('snps', snps)
    lines = []
    for s in range(len(snps)):
        lines.append([alt[snps[s]], str(int(line[1])+snps[s])])

    return lines

def process_dels(line):
    print('dels', line)
    #lines = []
    ref = line[3]
    alt = line[4]
    end = len(ref)
    assert len(alt) == 1
    if alt[0] == ref[0]:
        print(ref, alt, end)
        #gap or N? 
        #make sure the alignemnt is correct for this 
        l = ['-', str(int(line[1])+1), str(end-1)]

    
    elif alt[0]  == 'N':
        l = ['N', line[1], str(end)]

    else:
        l = line
        print('wtf', line)

    #print('lines', line)
    return l

#currently I will mask the whole region and then figure it out later 

def process_others(line):
    print('line', line)
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['N', line[1], str(end), '*']
    if ref > alt:
        l = ['N', line[1], str(end), '*']

    #print('weird',line)
    return l





        
        
    #create a new line that can be modified wo messing up line
    '''
    newLine = copy.copy(line)
    #remove allele count to avoid confusion (modify later if worth time)
    newLine[7] = 'AC='+';'+line[7][1]
    #iterate through all samples to update genotypes 
    for j in range(9, len(newLine)):
        #if genotype matches your current alt, make it 1
        if newLine[j] == str(a+1):
            newLine[j] = '1'
        else:
            newLine[j] = '0'
    #iterate through all snps and update newLine
    #add newLine to addLines before updating again 
    for snp in snps:
        newLine[1] = str(snp+ int(line[1]))
        newLine[3] = ref[snp]
        newLine[4] = newAlts[a][snp]
        #print('new',newLine)
        #yield line
        addLines.append(newLine)
    #if snp is the only alt, dont add line to newLines
    if len(newAlts) == 1:
        dontAdd = True
    else:
        removeAlt.append(a)
        '''
                    


def vcf_to_diff(vcf_file, output):

    with open(vcf_file) as v:
        with open(output, 'w') as o:
            missing = 0
            total = 0
            for line in v:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.strip().split()
                        sample = line[-1]
                        #print(sample)
                        o.write(f'>{sample}\n')
                    else:
                        total += 1
                        line = line.strip().split()
                        var = line[-1]
                        
                        if var != '0/0':
                            #print(line)
                            if var == './.':
                                print('missing', line)
                                missing += 1
                                line[4] = 'N'
                                line [-1] = '1'

                            else: 
                                print(line)
                                var = var.split('/')
                                var = var[0]
                                #print(var)
                                alts = line[4].split(',')
                                #if len(alts) > 1:
                                #    print(line[1],line[3],alts,var)
                                alt = alts[int(var)-1]
                                #print(alt)
                                line[4] = alt
                                line[-1] = '1'

                            

                            
                            #print('line', line)
                            assert type(line[4]) == str
                            if len(line[3]) == 1:

                                if len(line[4]) == 1:
                                    #print(line)
                                    o.write('\t'.join([line[4],line[1]])+'\n')
                                #elif len(line[4]) > 1:
                                    #print('insertion')


                            elif len(line[3]) > 1:
                                if len(line[4]) == len(line[3]):
                                    print(line)
                                    newlines = find_snps(line)
                                    for n in newlines:
                                        o.write('\t'.join(n)+'\n')

                                elif len(line[4]) == 1:
                                    #print('deletion')
                                    newline = process_dels(line)
                                    print(newline)
                                    o.write('\t'.join(newline)+'\n')

                                else:
                                    newline = process_others(line)
                                    o.write('\t'.join(newline)+'\n')



                                    #print('indel', line)
                                #print(line)
            print('missing', missing, 'total', total)

'''

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
                #print('lenRow',lenRow)
                break
            '''


#note this assumes 9 meta cols, will need to change here as well 
#this is all 1-indexing dont get confused




os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}filtvcf.vcf {vcf}")    
vcf_to_diff(f'{wd}filtvcf.vcf', 'output.vcf')