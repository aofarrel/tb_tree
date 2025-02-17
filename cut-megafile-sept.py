import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--megaVCF', required=True, type=str,help='path to giant VCF to be chunked')
parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-c', '--numMetadataColumns', default=9, type=int, help='one-indexed number of metadata columns in VCF, excluding sample names '
    '(ex: if CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SRR7592347 then enter 9)'
    'Note that ALT is always considered to be number 5.')
#parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)

args = parser.parse_args()
vcf = args.megaVCF
wd = args.workingDirectory
metacolumns = args.numMetadataColumns
#t = args.threads

#makes sure input path wont cause error
if wd[-1] != '/':
    #print('add a final slash')
    wd = wd+'/'
#sd = args.scriptsDirectory
#if sd[-1] != '/':
    #print('add a final slash')
#    sd = sd+'/'

                            
def find_snps(line):
    #for lines where len(ref)==len(alt), look for snps instead of processing as one large chunk

    ref = line[3]
    alt = line[4]
    end = len(ref)
    snps = []
    for i in range(end):
        if ref[i] != alt[i]:
            snps.append(i)
    #generate new lines for diff file
    lines = []
    for s in range(len(snps)):
        lines.append([alt[snps[s]], str(int(line[1])+snps[s]), '1'])

    return lines

def process_dels(line):
    #for lines where len(ref) > 1 AND len(alt) == 1, these deletions are easily processed 
    #currently deletions and missing data are all converted to '-'
    ref = line[3]
    alt = line[4]
    end = len(ref)
    assert len(alt) == 1

    #make sure remaining alt nucleotide is the same as the corresponding ref nuc
    if alt[0] == ref[0]:
        l = ['-', str(int(line[1])+1), str(end-1)]

    
    elif alt[0]  == '-':
        #this is for missing data
        l = ['-', line[1], str(end)]

    else:
        #this is a scenario that could be represented by a snp at the first ref position
        #this can be changed 
        l = ['-', line[1], str(len(line[3]))]
        
    return l



def process_others(line):
    #for sitations where reference and alt do not align 
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    return l


def vcf_to_diff(vcf_file, output):
    #takes a single sample vcf and converts to diff format 
    with gzip.open(vcf_file, 'rt') as v:
        with open(output, 'w') as o:
            missing = 0
            total = 0
            for line in v:
                if not line.startswith('##'):
                    if line.startswith('#'):
                        line = line.strip().split()
                        sample = line[-1]
                        #print('sample', sample)
                        o.write(f'>{sample}\n')
                    else:
                        total += 1
                        line = line.strip().split()
                        var = line[-1]
                        
                        if var != '0/0':
                            
                            if var == './.':
                                #print(f'missing {line}')
                                missing += 1
                                line[4] = '-'
                                line [-1] = '1'

                            else: 

                                var = var.split('/')
                                var = var[0]
                                alts = line[4].split(',')
                                alt = alts[int(var)-1]
                                line[4] = alt
                                line[-1] = '1'

                            #print(f'line {line}')
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
                                    #print(newline)
                                    o.write('\t'.join(newline)+'\n')

                                else:
                                    newline = process_others(line)
                                    o.write('\t'.join(newline)+'\n')



                                    #print('indel', line)
                                #print(line)
            print(f'missing {missing}, total {total}')

    return sample



def squish(file):
    #condenses diff entries that can be 
    with open(file) as f:
        with open(f'{file}squish', 'w') as o:

            prev = None
            for line in f:
                if not line.startswith('>'):

                    #process line
                    line = line.strip().split()
                    #print(line, prev)
                    
                    if prev != None:
                        #print('line', line)
                        #print('prev', prev)
                        if prev[0] == line[0]:
                            #print('prev', prev, 'line', line)
                            if int(prev[1]) == int(line[1])-int(prev[2]):
                                #print('prev', prev, 'now', line)
                                prev[2] = str(int(prev[2])+int(line[2]))
                                #print('new prev', prev)
                            else:
                                o.write('\t'.join(prev)+'\n') 
                                prev = line
                        else:
                            o.write('\t'.join(prev)+'\n') 
                            prev = line
                        
                    
                    #the first line of file becomes prev variable 
                    else:
                        prev = line
                        #print('first prev!', prev)
                        
                        
            

                #write header to new file
                else:
                    o.write(line)

                            
                            



#open file to indentify number of samples 
with gzip.open(vcf, 'rt') as v:
    for line in v:
        #if the line is not part of the heading
        if not line.startswith('##'):

            #if the line contains the column names 
            #note this assumes args.numMetadataColumns is set correctly
            if line.startswith('#'):
                line = line.strip().split('\t')
                #numSamps = len(line[metacolumns:])
                lenRow = len(line)

                #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                #print('lenRow',lenRow)
                break


#change this loop after testing 

assert(metacolumns != lenRow)
for i in range(metacolumns,lenRow):
    print(i)
    os.system(f'gzip -dc {vcf} | cut -f1-{lenRow},{i} > {wd}col{i}.vcf')
    os.system(f'bgzip -f {wd}col{i}.vcf')
    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {wd}col{i}filt.vcf {wd}col{i}.vcf.gz")
    os.system(f'bgzip -f {wd}col{i}filt.vcf')
    sample = vcf_to_diff(f'{wd}col{i}filt.vcf.gz', f'{wd}col{i}.diff')
    print(sample)
    os.system(f'rm {wd}col{i}*.vcf.gz')
    squish(f'{wd}col{i}.diff')
    os.system(f'rm {wd}col{i}.diff')
    if '/' not in sample:
        os.system(f'mv {wd}col{i}.diffsquish {wd}{sample}.diff')
    else:
        newname = sample.replace('/', '-')
        print(newname)
        os.system(f'mv {wd}col{i}.diffsquish {wd}{newname}.diff')
        

    #vcf_to_diff(f'{wd}col{i}filt.vcf', i, f'{wd}col{i}.diff')