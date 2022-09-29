import os
import argparse
import gzip
#import merged_to_diff_sept as mtd

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--megaVCF', required=True, type=str,help='path to giant VCF to be chunked')
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)
parser.add_argument('-sd', '--scriptsDirectory', type=str, required=True, help='path to directory where scripts are')

args = parser.parse_args()
vcf = args.megaVCF
wd = args.workingDirectory

#makes sure input path wont cause error
if wd[-1] != '/':
    #print('add a final slash')
    wd = wd+'/'
sd = args.scriptsDirectory
if sd[-1] != '/':
    #print('add a final slash')
    sd = sd+'/'

t = args.threads

                            
def find_snps(line):
    print('snps', line)
    ref = line[3]
    alt = line[4]
    end = len(ref)
    snps = []
    for i in range(end):
        if ref[i] != alt[i]:
            snps.append(i)
    
        
    #print('snps', snps)
    lines = []
    for s in range(len(snps)):
        lines.append([alt[snps[s]], str(int(line[1])+snps[s]), '1'])

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

    
    elif alt[0]  == '-':
        l = ['-', line[1], str(end)]

    else:
        l = ['-', line[1], str(len(line[3]))]
        print('wtf', line)

    #print('lines', line)
    return l

#currently I will mask the whole region and then figure it out later 

def process_others(line):
    print('others', line)
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    #print('weird',line)
    return l


def vcf_to_diff(vcf_file, output):

    with gzip.open(vcf_file, 'rt') as v:
        with open(output, 'w') as o:
            missing = 0
            total = 0
            for line in v:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.strip().split()
                        sample = line[-1]
                        print('sample', sample)
                        o.write(f'>{sample}\n')
                    else:
                        total += 1
                        line = line.strip().split()
                        var = line[-1]
                        
                        if var != '0/0':
                            
                            if var == './.':
                                print('missing', line)
                                missing += 1
                                line[4] = '-'
                                line [-1] = '1'

                            else: 
                                
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

                            

                            
                            print('line', line)
                            assert type(line[4]) == str
                            if len(line[3]) == 1:

                                if len(line[4]) == 1:
                                    #print(line)
                                    o.write('\t'.join([line[4],line[1], '1'])+'\n')
                                #elif len(line[4]) > 1:
                                    #print('insertion')


                            elif len(line[3]) > 1:
                                if len(line[4]) == len(line[3]):
                                    #print(line)
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

    return sample

                            
                            



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



#note this assumes 9 meta cols, will need to change here as well 
#this is all 1-indexing dont get confused

#change this loop after testing 

#for i in range(10,lenRow):
for i in range(10,16):
    print(i)

    os.system(f'gzip -dc {vcf} | cut -f1-9,{i} > {wd}col{i}.vcf')
    os.system(f'bgzip -f {wd}col{i}.vcf')
    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}col{i}filt.vcf {wd}col{i}.vcf.gz")
    os.system(f'bgzip -f {wd}col{i}filt.vcf')
    sample = vcf_to_diff(f'{wd}col{i}filt.vcf.gz', f'{wd}col{i}.diff')
    os.system(f'mv {wd}col{i}.diff {wd}{sample}.diff')
    #vcf_to_diff(f'{wd}col{i}filt.vcf', i, f'{wd}col{i}.diff')