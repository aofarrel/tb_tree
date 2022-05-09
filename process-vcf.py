import argparse
import gzip as gz

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcfFile', required=True)
parser.add_argument('-o', '--outFile', required=True)
args = parser.parse_args()
inFile=args.vcfFile
outFile=args.outFile

def processIndels(line, alleleCounts, alts, variantPresent, missingOnly, outFile):
    print(line)
    ref = line[3]
    print('ref', ref)
    print(alleleCounts, alts, variantPresent, missingOnly)
    newAlts = []
    newAltsRef = []
    newAC = []
    if missingOnly == False:
        for i in range(len(alleleCounts)):
            print(alleleCounts[i])
            if int(alleleCounts[i]) > 0:
                #make sure ref is 1 indexed
                newAltsRef.append(i+1)
                newAlts.append(alts[i])
                newAC.append(alleleCounts[i])
        print('new ac', newAC, 'newAltsRef', newAltsRef, 'newAlts', newAlts)
        assert len(newAltsRef) == len(newAlts)
        assert len(newAlts) == len(newAC)

        line[7].split(';')[0] = 'AC='+ ','.join(newAC)
        print('new',line)




    else:
        print('what are you?', line)

    #for i in range(9,len(line)):
    #    print(line[i])

#def processSNPs

with gz.open(inFile, 'rb') as f:
    with open(outFile, 'w') as of:
        deleted = 0
        missingData = 0
        detected = 0
        processed = 0
        lines = 0
        indels = 0 
        for line in f:
            
            #convert from bytes to str 
            line = str(line, 'UTF-8')
            #skip through header
            if line.startswith('#'):
                of.write(line)
                continue
            
            lines += 1
            #cut up line
            line = line.strip().split('\t')
            #note that alleleCounts are 2X the amount bc it assumes diploid
            alleleCounts = line[7].split(';')[0][3:].split(',')
            #print('AC',alleleCounts)
            alts = line[4].split(',')
            #print('alts',alts)
            #make sure there is a count for every allele
            assert len(alts) == len(alleleCounts)
            
            #iterate through alleleCounts to determine if there are variants at the position
            variantPresent = False
            for count in alleleCounts:
                if int(count) > 0:
                    variantPresent = True
                    break
            #print('variant? ', variantPresent)

            # if there are no variants, check for missing data before deleting 
            onlyMissing = False
            if variantPresent==False:
                if './.' not in line:
                    #there is no missing data or variants, this line can be ignored
                    deleted += 1
                    continue
                else:
                    #there is missing data in this line and it must be processed
                    missingData += 1
                    onlyMissing = True
            else:
                detected += 1
            missing = False
            missingperline = 0

            # identify if line contains indels 
            #indels only need to be processed if there are variants in the line (if there is only missing data they can be processed as SNPs)?
            #if onlyMissing == False:
            possibleIndel = False
            for a in alts:
                if len(a) > 1:
                    possibleIndel = True
                    break

            if possibleIndel == True or len(line[3]) > 1:
                print('indel?')
                processIndels(line, alleleCounts, alts, variantPresent, onlyMissing, of)
                indels += 1
            
            #print('indel', indel)

            # iterate through all samples at position to find missing data, change genotypes from diploid to haploid and write to file
            else:
                for i in range(9,len(line)):
                    gt = line[i].split('/')
                    #when written as diploid the gt should have matching haplotype
                    assert gt[0] == gt[1]
                    #if a sample has missing data, determine if 'N' needs to be added to alt alleles and change genotype to the index for 'N'
                    if line[i] == './.':
                        missingperline += 1
                        if missing == False:
                            line[4] += ',N'
                            missing = True
                        
                        #check this!!!!
                        #if ref is 0 then all alts are 1 based encoding?
                        line[i] = str(line[4].split(',').index('N')+1)
                        
                    #if no data is missing for sample, change genotype to haploid   
                    else:
                        line[i] = gt[0]

                #print('processed', line)
                #note that AC does not count the number of N alleles per line 
                #print('missingperline', missingperline)
                of.write('\t'.join(line)+'\n')
                processed += 1
                
        print('number of deleted positions', deleted)
        print('positions w missing data', missingData )
        print('variants detected', detected)
        print('positions processed (detected + missingData)', processed)
        print('indels', indels)
        print('total', lines)
