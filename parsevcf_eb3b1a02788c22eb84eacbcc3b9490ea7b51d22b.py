# see https://github.com/lilymaryam/parsevcf/commit/eb3b1a02788c22eb84eacbcc3b9490ea7b51d22b

import os
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--VCF', required=True, type=str,help='path to VCF to be processed')
parser.add_argument('-d', '--working_directory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-tbmf', '--tb_maskfile', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-cf', '--bed_coverage_file', required=False, type=str, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")
parser.add_argument('-cd', '--coverage_depth', required=False, default=10, type=int, help="path to bed coverage file for vcf (note: can only be used with single-sample vcfs)")

args = parser.parse_args()
vcf = args.VCF
wd = args.working_directory
tbmf = args.tb_maskfile
cf = args.bed_coverage_file
cd = args.coverage_depth

#makes sure input path wont cause error
if wd[-1] != '/':
    wd = wd+'/'

#Functions

#for lines where len(ref)==len(alt), look for snps instead of processing as one large chunk                        
def find_snps(line):
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

#for lines where len(ref) > 1 AND len(alt) == 1 
#currently deletions and missing data are all converted to '-'
def process_dels(line):
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
        l = ['-', line[1], str(len(line[3]))]
        
    return l


#for sitations where reference and alt do not align
#will mask entire reference 
def process_others(line):
    ref = line[3]
    alt = line[4]
    end = len(ref)
    if ref < alt:
        l = ['-', line[1], str(end)]
    if ref > alt:
        l = ['-', line[1], str(end)]

    return l

#read TB mask file (different than coverage file) and generate sites to be masked
#bed coverage file is 0 indexed, so add 1 to everything 
def mask_TB(tbmf):
    tb_sites = {}
    with open(tbmf) as file:
        for line in file:
            line=line.strip().split()
            tb_sites[int(line[1])+1] = int(line[2])+1
    #print('regions', tb_sites)
    return tb_sites

#read coverage file and generate sites to be masked 
#if coverage does not have HR37c reference it will throw an error (this can be changed)
#bed files are 0-indexed in col1 and 1-indexed in col2, i am adding one to both to make them one indexed
#and to maintain their mapping logic
def mask_low_depth(cf, cd):
    ld_sites = {}
    prev = None
    with open(cf) as cf:
        for line in cf:
            #for currect bed coverage file 
            if line.startswith('NC_000962.3'):
                line = line.strip().split()
                #print('pre',line)
                line[1] = str(int(line[1])+1)
                line[2] = str(int(line[2])+1)
                #print('post',line)
                if int(line[3]) < cd:
                    #print('prev', prev)
                    #print(line)
                    if prev == None:
                        ld_sites[int(line[1])] = int(line[2])
                        prev = [int(line[1]), int(line[2])]
                    else:
                        if int(line[1]) == prev[1]:
                            ld_sites[prev[0]] = int(line[2])
                            #print('squish', 'prev', prev, 'line', line)
                            prev[1] = int(line[2])

                        else:
                            #print('no squish')
                            ld_sites[int(line[1])] = int(line[2])
                            prev = [int(line[1]), int(line[2])]
                    

    if ld_sites == {}:
        raise Exception('coverage file has incorrect reference')
    #print('regions', count)
    #print('ld_sites', len(ld_sites))
    #for l in ld_sites:
    #    print(l, ld_sites[l])
    return ld_sites

#when merging ld and tb masks need to make sure the added masks are not overlapping
def check_prev_mask(prev, line):
    #currently not checking overlap to left of prev bc that indicates a bigger error
    #may need to change?
    overlap = False
    change = None
    #print('prev', prev)
    #print('line', line)
    #if prev != None:
    prev_s = int(prev[0])
    prev_e = int(prev[1])
    line_s = int(line[0])
    line_e = int(line[1])
    #if prev != None:
    #print('prev', prev_s, prev_e, 'line', line_s, line_e)
    if line_s >= prev_s and line_e <= prev_e:
        #print('prev', prev_s, prev_e, 'line', line_s, line_e)
        overlap = True
        #print('Full OVERLAP!!!!!')
        #print('overlap',overlap, 'change', change)
    elif line_s >= prev_s and line_s <= prev_e and line_e >= prev_e:
        #print('prev', prev_s, prev_e, 'line', line_s, line_e)
        overlap = True 
        #print('right overlap!!!!!!')
        prev[1] = line_e
        change = prev
        #print('overlap',overlap, 'change', change)
    return overlap, change


def condense_mask_regions(cf,cd,tbmf):
    ld_sites = mask_low_depth(cf, cd)
    tb_sites = mask_TB(tbmf)
    tb_keys = sorted(tb_sites.keys())
    ld_keys = sorted(ld_sites.keys())
    all_sites = {}
    #print('ld',ld_keys)
    #print('tb',tb_keys)
    tb_keys_ind = 0
    ld_keys_ind = 0
    cont = 0
    prev = None
    while tb_keys_ind < len(tb_keys) or ld_keys_ind < len(ld_keys):
        
            #print(all_sites_keys[-1])
        if tb_keys_ind < len(tb_keys) and ld_keys_ind < len(ld_keys): 
            tb_start = tb_keys[tb_keys_ind]
            tb_end =  tb_sites[tb_keys[tb_keys_ind]]
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]

            if ld_start < tb_start:
                if ld_end < tb_start:
                    #print(f'ld{ld_keys_ind} is below tb{tb_keys_ind}')
                    #print('ld', ld_start, ld_end)
                    #print('tb', tb_start, tb_end)
                    if prev != None:
                        overlap, change = check_prev_mask(prev, [ld_start, ld_end])
                        if overlap == True and change != None:
                            #print('change', 'need to update prev!!!!!', change)
                            #print(f'ld{ld_keys_ind} is below tb{tb_keys_ind}')
                            #print('ld', ld_start, ld_end)
                            #print('tb', tb_start, tb_end)
                            all_sites[change[0]] = change[1]
                    if prev == None or overlap == False:
                        all_sites[ld_start] = ld_end
                    ld_keys_ind += 1
                elif ld_end >= tb_start:
                    #print('left overlap')
                    #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                    if prev != None:
                        overlap, change = check_prev_mask(prev, [ld_start, tb_end])
                        if overlap == True and change != None:
                            #print('prev', prev, 'change', change)
                            all_sites[change[0]] = change[1]
                    if prev == None or overlap == False:
                        #if prev == None:
                        #    print('first')
                        #elif overlap == False:
                        #    print('overlap = ', overlap)
                        all_sites[ld_start] = tb_end
                    tb_keys_ind += 1
                    ld_keys_ind += 1
            
            elif ld_start <= tb_end and ld_end > tb_end:
                #print('right over lap')
                #print('left overlap')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                if prev != None:
                    overlap,change = check_prev_mask(prev, [tb_start, ld_end])
                    if overlap == True and change != None:
                        #print('change', change)
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[tb_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            elif ld_start > tb_end:
                #print('no overlap')
                #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)
                if prev != None:
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != None:
                        #print('change', change)
                        #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                        #print('ld', ld_start, ld_end)
                        #print('tb', tb_start, tb_end)
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1

            elif ld_start >= tb_start and ld_end <= tb_end:
                #print('full overlap: ld inside')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)
                if prev != None:
                    overlap,change = check_prev_mask(prev, [tb_start, tb_end])
                    if overlap == True and change != False:
                        #print('change', change)
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[tb_start] = tb_end
                tb_keys_ind += 1
                ld_keys_ind += 1

            
            elif ld_start <= tb_start and ld_end >= tb_end:
                #print('full overlap: tb inside')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                if prev != None:
                    overlap,change = check_prev_mask(prev, [ld_start, ld_end])
                    if overlap == True and change != None:
                        #print('change', change)
                        all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    all_sites[ld_start] = ld_end
                tb_keys_ind += 1
                ld_keys_ind += 1
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)

            #else:
                #print('other', 'what else could happen?')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                #all_sites[tb_start]


                #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)
            #print('tb',tb_keys[tb_keys_ind], tb_sites[tb_keys[tb_keys_ind]])
            #print('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
        elif tb_keys_ind >= len(tb_keys) and ld_keys_ind < len(ld_keys):
            #print('no more tb masks, ld only')
            ld_start = ld_keys[ld_keys_ind]
            ld_end = ld_sites[ld_keys[ld_keys_ind]]
            all_sites[ld_start] = ld_end
            #tb_keys_ind += 1
            ld_keys_ind += 1


        elif ld_keys_ind >= len(ld_keys) and tb_keys_ind < len(tb_keys):
            #print('no more ld masks, tb only')
            tb_start = tb_keys[tb_keys_ind]
            tb_end = tb_sites[tb_keys[tb_keys_ind]]
            all_sites[tb_start] = tb_end
            #tb_keys_ind += 1
            tb_keys_ind += 1 

        cont += 1
        #print('tb ind', tb_keys_ind)
        #print('ld ind', ld_keys_ind)
        #print('len tb', len(tb_keys))
        #print('len ld', len(ld_keys))
        #if cont == 10000:
        #    print(all_sites)
        #    break
        #prev = 
        if len(all_sites) > 0:
            all_sites_keys = sorted(all_sites.keys())
        prev = [all_sites_keys[-1],all_sites[all_sites_keys[-1]]]
    return all_sites
                    
def squish(lines):
    #condenses diff lines that can be compressed into a single line
    #print('SQUISH')
    prev = None
    newLines = []
    for line in lines:
        #print(line)
        if not line[0].startswith('>'):
            if prev != None:
                #print('line', line)
                #print('prev', prev)
                if prev[0] == line[0]:
                    
                    if int(prev[1]) == int(line[1])-int(prev[2]):
                        #print('squish', 'prev', prev, 'line', line)
                        #print('prev', prev, 'now', line)
                        prev[2] = str(int(prev[2])+int(line[2]))
                        #print('new prev', prev)
                    else:
                        #o.write('\t'.join(prev)+'\n')
                        newLines.append(prev) 
                        prev = line
                else:
                    #o.write('\t'.join(prev)+'\n')
                    newLines.append(prev) 
                    prev = line
                    
                
            
            #the first line of file becomes prev variable 
            else:
                
                #print('first prev!', prev)
                #newLines.append(line)
                prev = line
        #write header to new file
        else:
            newLines.append(line)
            #print(line)
            #o.write(line) 
    if prev != None:
        newLines.append(prev)  
    return newLines

                              
#option to write output to file
def vcf_to_diff(vcf_file, output):
    #takes a single sample vcf and converts to diff format 
    lines = []
    with open(vcf_file, 'rt') as v:
        #with open(output, 'w') as o:
        missing = 0
        total = 0
        for line in v:
            if not line.startswith('##'):
                if line.startswith('#'):
                    line = line.strip().split()
                    sample = line[-1]
                    #print('sample', sample)
                    #o.write(f'>{sample}\n')
                    lines.append([f'>{sample}'])
                else:
                    total += int(len(line[3]))
                    line = line.strip().split()
                    var = line[-1]
                    
                    if var != '0/0':
                        
                        if var == './.':
                            #print('missing', line)
                            #missing += int(len(line[3]))
                            line[4] = '-'
                            line [-1] = '1'

                        else: 
                            
                            var = var.split('/')
                            var = var[0]
                            alts = line[4].split(',')
                            alt = alts[int(var)-1]
                            line[4] = alt
                            line[-1] = '1'

                        #print('line', line)
                        assert type(line[4]) == str
                        if len(line[3]) == 1:

                            if len(line[4]) == 1:
                                #print(line)
                                #o.write('\t'.join([line[4],line[1], '1'])+'\n')
                                lines.append([line[4],line[1], '1'])
                            #elif len(line[4]) > 1:
                                #print('insertion')


                        elif len(line[3]) > 1:
                            if len(line[4]) == len(line[3]):
                                #print(line)
                                newlines = find_snps(line)
                                for n in newlines:
                                    #o.write('\t'.join(n)+'\n')
                                    lines.append(n)

                            elif len(line[4]) == 1:
                                #print('deletion')
                                newline = process_dels(line)
                                #print(newline)
                                #o.write('\t'.join(newline)+'\n')
                                lines.append(newline)

                            else:
                                newline = process_others(line)
                                #o.write('\t'.join(newline)+'\n')
                                lines.append(newline)
    #figure out how to count missing data and add stats here 
    #with open('testlines','w') as o:
    #    for line in lines:
    #        o.write('\t'.join(line)+'\n')
    newLines = squish(lines)
    return newLines

                                    #print('indel', line)
                                #print(line)
            #print('missing', missing, 'total', total)

    #return sample, missing, total
    #return sample, missing, total


def make_files(lenRow, samps,wd):
    files = {}
    for i in range(len(samps)):
        #print(i)
        s = samps[i]
        #print(s)
        if '/' not in s:
            file = open(f'{wd}{s}.vcf','w')
        else:
            newname = s.replace('/', '-')
            #print(newname)
            file = open(f'{wd}{newname}.vcf','w')
        #files.append(file)
        files[i+9] = file
        #print(files)
    return files

                            
def count_samples(vcf):
    #open file to indentify number of samples 
    with open(vcf, 'r') as v:
        for line in v:
            #if the line is not part of the heading
            if not line.startswith('##'):

                #if the line contains the column names 
                #note that this assumes 9 meta data columns, if you need to account for more/less change this 
                if line.startswith('#'):
                    line = line.strip().split('\t')
                    samps = line[9:]
                    lenRow = len(line)

                    #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                    # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                    #print('lenRow',lenRow)
                    break
    
    return lenRow, samps                           


def count_samples_bin(vcf):
    #open file to indentify number of samples 
    with gzip.open(vcf, 'rt') as v:
        for line in v:
            #if the line is not part of the heading
            if not line.startswith('##'):

                #if the line contains the column names 
                #note that this assumes 9 meta data columns, if you need to account for more/less change this 
                if line.startswith('#'):
                    line = line.strip().split('\t')
                    samps = line[9:]
                    lenRow = len(line)

                    #this is the 1-indexed indices for each row of the file, we will keep these values as they can be 
                    # traced from the largest vcf to the subvcfs to make sure all of the info is consistent
                    #print('lenRow',lenRow)
                    break
    
    return lenRow, samps

def read_VCF(vcf, files):
    #count = 0 
    with open(vcf) as v:
        for line in v:
            if not line.startswith('##'):
                line = line.strip().split()
                position = line[0:9] 
                #print('position', position)
                #print(len(line))
                #count = 1
                for f in files:
                    #print('f',f)
                    #print(files[f])
                    #count += 1
                    parcel = [line[f]]
                    #print('parcel',parcel)
                    #print(parcel)
                    newline = position+parcel
                    #print(newline)
                    #print(len(newline))
                    files[f].write('\t'.join(newline)+'\n')
                    #print(position[1],f,files[f][0])
                    #print(f)
                    
                #count += 1

            
                
            
            else:
                #print(line)
                for f in files:
                    files[f].write(line)

def read_VCF_bin(vcf, files):
    #count = 0 
    with gzip.open(vcf, 'rt') as v:
        for line in v:
            if not line.startswith('##'):
                line = line.strip().split()
                position = line[0:9] 
                #print('position', position)
                #print(len(line))
                #count = 1
                for f in files:
                    #print('f',f)
                    #print(files[f])
                    #count += 1
                    parcel = [line[f]]
                    #print('parcel',parcel)
                    #print(parcel)
                    newline = position+parcel
                    #print(newline)
                    #print(len(newline))
                    files[f].write('\t'.join(newline)+'\n')
                    #print(position[1],f,files[f][0])
                    #print(f)
                    
                #count += 1

            
                
            
            else:
                #print(line)
                for f in files:
                    files[f].write(line)

    #else:
    #    print('prev', prev, 'line', line_s, line_e)

def check_prev_line(prev, line):
    #currently not checking overlap to left of prev bc that indicates a bigger error
    #may need to change?
    overlap = False
    change = None
    print('prev', prev)
    print('line', line)
    #if prev != None:
    prev_s = int(prev[1])
    prev_e = prev_s + int(prev[2])
    line_s = int(line[1])
    line_e = line_s + int(line[2])
    #if prev != None:
    
    print('prev', prev_s, prev_e, 'line', line_s, line_e)
    if line_s >= prev_s and line_e <= prev_e:
        overlap = True
        print('Full OVERLAP!!!!!')
        print('prev', prev_s, prev_e, 'line', line_s, line_e)
    elif line_s >= prev_s and line_s <= prev_e and line_e >= prev_e:
        overlap = True 
        print('right overlap!!!!!!')
        print('prev', prev)
        print('line', line)
        print('prev', prev_s, prev_e, 'line', line_s, line_e, line)
        if line_s == prev_e and line[-1]=='1':
            #print('SNP')
            overlap = False
            #change = line
        elif line_s == prev_e and line[-1]!='1':
            overlap = True
            print('need to change this!!!', 'prev', prev, 'line', line)
            print('prev')
        else:
            print('what is happening?')
            print(f'{line_e}-{prev_s}=',line_e-prev_s)
            print('prev', prev)
            prev[2] = str(int(line[2])-int(prev[1]))
            print('prev after ', prev)
            change = prev
    #elif line_s <= prev_s and line_e >= prev_s:
    print('overlap',overlap, 'change', change)
    return overlap, change
    

def mask_and_write_diff(cf, cd, tbmf, lines, samps, sample):
    #print(cf)
    #print(cd)
    #print(tbmf)
    #print(lines)
    #print(samps)
    #print('sample', sample)
    
    if cf != None:
        assert len(samps) == 1
        masks = condense_mask_regions(cf,cd,tbmf)
        with open('masks','w') as f:

            for m in masks:
                f.write('\t'.join([str(m),str(masks[m])])+'\n')
        
    else:
        #print('no coverage mask file')
        masks = mask_TB(tbmf)

    #print('masks',len(masks))
    #print('lines', len(lines))
    masks_key = sorted(masks.keys())
    #ld_keys = sorted(ld_sites.keys())
    masks_ind = 0
    lines_ind = 1


    '''
    prev = None
    for line in lines:
        if len(line) == 1:
            print(line)
        else:
            print('prev', prev)
            if prev != None:
                if int(line[1]) <= prev:
                    print('OUT OF ORDER ################################')
                else:
                    print('prev', prev)
                    print('line',line)
                    #pass
            
            prev = int(line[1])
            '''
    all_lines = [lines[0]]
    cont = 0
    prev = None
    while masks_ind < len(masks_key) or lines_ind < len(lines):
        #print('prv', prev)
        #if len(all_lines) > 1:
        #    all_lines_keys = sorted(all_lines.keys())
        
        #if both indexes are still going
        if masks_ind < len(masks_key) and lines_ind < len(lines): 
            mask_start = masks_key[masks_ind]
            mask_end =  masks[mask_start]
            line = lines[lines_ind]
            line_start = int(lines[lines_ind][1])
            line_end = int(lines[lines_ind][1])+int(lines[lines_ind][2])
            print('mask_start', mask_start, type(mask_start))
            print('mask_end', mask_end, type(mask_end))
            print('line',line)
            print('line_start', line_start, type(line_start))
            print('line end', line_end, type(line_end))

            if line_start >= mask_start and line_end <= mask_end:
                print('full overlap: line inside')
                print('line',line_start, line_end, 'mask', mask_start, mask_end)
                print('line', line_start, line_end)
                print('mask', mask_start, mask_end)
                

                if prev != None:
                    overlap,change = check_prev_line(prev, ['-', mask_start, mask_end])
                    if overlap == True and change != None:
                        print('change', change)
                        #all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    
                    all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])


                #if prev != None:
                #    check_prev_line(prev, line)
                #print('prev', prev)
                #print('append?', ['-', str(mask_start), str(mask_end-mask_start)])
                #
                #all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
                masks_ind += 1
                lines_ind += 1

            #need to make sure that if snps overlap they get masked
            elif line_start <= mask_start and line_end >= mask_end:
                print('full overlap: mask inside')
                print('line',line_start, line_end, 'mask', mask_start, mask_end)
                print('prev', prev)
                print('append?', line)

                if prev != None:
                    overlap,change = check_prev_line(prev, line)
                    if overlap == True and change != None:
                        print('change', change)
                        #all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                    
                    all_lines.append(line)


                #if prev != None:
                #    check_prev_line(prev, line)
                #if prev == None:
                #all_lines.append(line)
                masks_ind += 1
                lines_ind += 1
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)

        
            elif line_start < mask_start:
                #no overlap, line completely to left
                if line_end < mask_start:
                    print(f'ld{lines_ind} is below tb{masks_ind}')
                    print('no overlap')
                    #print('ld', ld_start, ld_end)
                    #print('tb', tb_start, tb_end)
                    print('prev', prev)
                    print('append?', line)

                    if prev != None:
                        overlap,change = check_prev_line(prev, line)
                        if overlap == True and change != None:
                            print('change', change)
                            print('before',all_lines[-1])
                            all_lines[-1][2] = change[1]
                            print('after',all_lines[-1])
                    if prev == None or overlap == False:
                    
                        all_lines.append(line)


                    #if prev != None:
                    #    check_prev_line(prev, line)
                    #if prev == None:
                    #all_lines.append(line)
                    lines_ind += 1
                
                elif line_end >= mask_start:
                    print('left overlap ##########################################')
                    print('line',line_start, line_end, 'tb', mask_start, mask_end)
                    print(line)
                    print('prev', prev)
                    #print('append?', ['-', str(mask_start), str(line_end)])


                    if prev != None:
                        overlap,change = check_prev_line(prev, line)
                        if overlap == True and change != None:
                            print('change', change)
                            #all_sites[change[0]] = change[1]
                    if prev == None or overlap == False:
                    
                        all_lines.append(['-', str(line_start), str(mask_end-line_end)])


                    #if prev != None:
                    #    check_prev_line(prev, line)
                    #if prev == None:
                    #all_lines.append(['-', str(mask_start), str(line_end)])
                    masks_ind += 1
                    lines_ind += 1
                    #all_ = tb_end
                    #tb_keys_ind += 1
                    #ld_keys_ind += 1

            elif line_start > mask_end:
                print('no overlap')
                #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                print('line', line_start, line_end)
                print('mask', mask_start, mask_end)
                print(f'["-", {str(mask_start)}, {str(mask_end-mask_start)}]')
                print('prev', prev)
                print('append?', ['-', str(mask_start), str(mask_end-mask_start)])

                if prev != None:
                    overlap,change = check_prev_line(prev, ['-',str(mask_start),str(mask_end)])
                    if overlap == True and change != None:
                        print('change', change)
                        #all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                
                    all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])


                #if prev != None:
                #    check_prev_line(prev, line)
                #if prev == None:
                #all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
                masks_ind += 1
            
            #need to figure out whatn happens if snp is sticking out 
            elif line_start <= mask_end and line_end > mask_end:
                print('right over lap ########################################################')
                print('ld',line_start, line_end, 'tb', mask_start, mask_end)
                print('prev', prev)
                #print('append?', ['-', str(mask_start), str(line_end-mask_start)])

                if prev != None:
                    overlap,change = check_prev_line(prev, line)
                    if overlap == True and change != None:
                        print('change', change)
                        #all_sites[change[0]] = change[1]
                if prev == None or overlap == False:
                
                    all_lines.append(['-', str(mask_start), str(line_end-mask_start)])


                #if prev != None:
                #    check_prev_line(prev, line)
                #if prev == None:
                #all_lines.append(['-', str(mask_start), str(line_end-mask_start)])
                masks_ind += 1
                lines_ind += 1
            
            
            
            
            else:
                print('other', 'what else could happen?')
                #print('ld',ld_start, ld_end, 'tb', tb_start, tb_end)
                #all_sites[tb_start]


                #print(f'ld{ld_keys_ind} is not below tb{tb_keys_ind}')
                #print('ld', ld_start, ld_end)
                #print('tb', tb_start, tb_end)
            #print('tb',tb_keys[tb_keys_ind], tb_sites[tb_keys[tb_keys_ind]])
            #print('ld',ld_keys[ld_keys_ind], ld_sites[ld_keys[ld_keys_ind]])
        
        elif masks_ind >= len(masks_key) and lines_ind < len(lines):
            print('no more tb masks, ld only')
            line = lines[lines_ind]
            line_start = line[1]
            #ld_end = int(line[2])+int(line[1])
            all_lines.append(line)
            #tb_keys_ind += 1
            lines_ind += 1


        elif lines_ind >= len(lines) and masks_ind < len(masks_key):
            #raise Exception('more masks ') 
            #print('no more ld masks, tb only')
            mask_start = masks_key[masks_ind]
            mask_end =  masks[mask_start]
            #mask_start = [tb_keys_ind]
            #tb_end = tb_sites[tb_keys[tb_keys_ind]]
            all_lines.append(['-', str(mask_start), str(mask_end-mask_start)])
            #tb_keys_ind += 1
            masks_ind += 1 
        
        cont += 1
        #print('tb ind', tb_keys_ind)
        #print('ld ind', ld_keys_ind)
        #print('len tb', len(tb_keys))
        #print('len ld', len(ld_keys))
        #print(cont)
        prev = all_lines[-1]
        print('prev', prev)
        #if cont == 100000:
        #    print(all_lines)
        #    break

    #prev = None
    #for l in all_lines:
    #    print(l)
        '''
        if prev == None:
            print('first')
        else:
            if int(l[1]) <= int(prev[2]):
                print('????', prev, l)
        
        prev = line
        '''
    #print('all lines', all_lines)
    #with open('diff.diff','w') as d:
    #    for line in all_lines:
    #        d.write('\t'.join(line)+'\n')
    return all_lines
    




#SCRIPT STARTS HERE

binary = True
with gzip.open(vcf, 'r') as test:
    try:
        test.read(1)
    except OSError:
        #print('not binary')
        binary = False

if binary == True:
    lenRow, samps = count_samples_bin(vcf)
else:
    lenRow, samps = count_samples(vcf)

#be careful w dictionaries!!!
files = make_files(lenRow, samps, wd)

if binary == True:
    read_VCF_bin(vcf, files)
    
else:
    read_VCF(vcf, files)
    
masks = mask_TB(tbmf)
#this is not parallelized, the more samples in the vcf the longer this will take
#print(files)
#print(masks)

for f in files:
    files[f].close()
    #print(f, files[f])
    
    sample = os.path.basename(files[f].name)[:-4]
    filepath = files[f].name
    #print('sample', sample)
    #print('filepath', filepath)

    
    #note: theres no point to do this since the files are already completely filling the space
    #    os.system(f'bgzip -f {filepath}')
    #    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}.gz")
    #    os.system(f"rm {filepath}.gz")
    
    
    os.system(f"bcftools annotate -x '^FORMAT/GT' -O v -o {filepath}.filt {filepath}")
    os.system(f"rm {filepath}")
    
    #note: theres no point to do this since the files are already completely filling the space
    #os.system(f'bgzip -f {filepath}.filt')
    #vcf_to_diff(f'{filepath}.filt.gz', f'{wd}{sample}.diff')

    #if there is a provided coverage file it will be used to mask low coverage (less than cd) regions 
    #note that only one coverage file can be provided and it will result in an error if the vcf has more samples than coverage files 

    #more work on this later 
    #print('FILTERING UNIVERSAL SITES ONLY')
    lines = vcf_to_diff(f'{filepath}.filt', f'{wd}{sample}.diff')
    os.system(f'rm {filepath}.filt')
    all_lines = mask_and_write_diff(cf,cd,tbmf,lines, samps, sample)
    with open(f'{wd}{sample}.diff','w') as o:
        for line in all_lines:
            #print(line)
            o.write('\t'.join(line)+'\n')
    #print(lines)
        #squish(f'{wd}{sample}.diff')
        #os.system(f'mv {wd}{sample}.diffsquish {wd}{sample}.diff')
    
    #figure out how to count missing samples, mask full low coverage regions

    #take out low depth masks see if it works 
    #else:
    #    print('FILTERING FOR LOW COVERAGE')
    #    assert len(samps) == 1, f'must only have one sample for each coverage file'
    #    lines = vcf_to_diff(f'{filepath}.filt', f'{wd}{sample}.diff')
    #    os.system(f'rm {filepath}.filt')
        #squish(f'{wd}{sample}.diff')
        #os.system(f'mv {wd}{sample}.diffsquish {wd}{sample}.diff')
        #all_masks = condense_mask_regions(cf,cd,tbmf)
        #for m in all_masks:
        #    print(m, all_masks[m])
        #print(len(all_masks))
        #mask_LDsites(cf,cd,f'{wd}{sample}.diff',sample)
        #os.system(f'mv {wd}{sample}.masked.diffsquish {wd}{sample}.diff')
    

    #os.system(f"rm {filepath}.filt")
    #mask_TBsites(f'{wd}{sample}.diff', tbmf, sample)
    #squish(f'{wd}{sample}.masked.diff')
    #os.system(f'mv {wd}{sample}.masked.diffsquish {wd}{sample}.masked.diff')
    #os.system(f'rm {wd}{sample}.diff')
    

#SCRIPT ENDS HERE





#ld = mask_low_depth(cf,cd)


'''
with open('testnewlines','w') as o:
    for line in lines:
        #print(line)
        o.write('\t'.join(line)+'\n')
        #print(line)
        
'''   
#mask_low_depth(cf,cd)
#all_masks = condense_mask_regions(cf,cd,tbmf)
#print(len(all_masks))

#for m in all_masks:
#    print(m,all_masks[m])


'''
prev = None
for m in all_masks:
    #print(m,all_masks[m])
    if prev == None:
        print('first')
    else:
        if m <= prev[1]:
            print('????', prev, m, all_masks[m])
    
    prev = [m, all_masks[m]]
    '''
    
    

    

    







