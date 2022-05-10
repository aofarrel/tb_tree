import argparse
import gzip as gz

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcfFile', required=True)
parser.add_argument('-o', '--outFile', required=True)
args = parser.parse_args()
inFile=args.vcfFile
outFile=args.outFile

def processAlleles(line):
    line = line.strip().split('\t')
    
    #note that alleleCounts are 2X the amount bc it assumes diploid
    alleleCounts = line[7].split(';')[0][3:].split(',')
    alts = line[4].split(',')

    #make sure there is a count for every allele
    assert len(alts) == len(alleleCounts)
    variantsPresent = False
    newAlts = []
    newAltsRef = []
    newAC = []
    for i in range(len(alleleCounts)):
        if int(alleleCounts[i]) > 0:
            #make sure ref is 1 indexed
            newAltsRef.append(str(i+1))
            newAlts.append(alts[i])
            newAC.append(alleleCounts[i])
            variantsPresent = True
    #make sure all new stats are equivalent 
    assert len(newAltsRef) == len(newAlts)
    assert len(newAlts) == len(newAC)

    onlyMissing = False
    missing = False
    variantOnly = False
    # if there are no variants, check for missing data before deleting 
    if variantsPresent==False:
        if './.' not in line:
            #there is no missing data or variants, this line can be ignored
            return None
        else:
            #there is missing data in this line and it must be processed
            onlyMissing = True
    
    else:
        if './.' in line:
            missing = True
        else:
            variantOnly = True
    
    assert onlyMissing == True or len(newAC) > 0
    
    #change allele counts (counts changed from dip to hap)
    #if greater efficiency is needed, dip->hap is not strictly necessary
    line[7] = line[7].split(';')
    changed = False
    one = 0
    two = 0
    three = 0

    #the first conditional accounts for all lines where onlyMissing = True
    if len(newAC) == 0:
        #print(line[1])
        newAlts.append('N')
        line[4] = 'N'
        line[7] = 'AC='+';'+line[7][1]
        #print('only missing', onlyMissing)
        #print('missing', missing)
        changed = True
        #print('1')
    elif alleleCounts == newAC:
        #print(line[1], 'unchanged')
        #print('missingOnly', onlyMissing)
        #print('missing', missing)
        #print('variantOnly', variantOnly)
        if missing == True:
            newAlts.append('N')
            line[4] = ','.join(newAlts)
        ac = line[7][0][3:].split(',')
        for a in range(len(ac)):
            ac[a] = str(int(int(ac[a])/2))
        line[7] = 'AC='+ ','.join(ac)+';'+line[7][1]
        #print('only missing', onlyMissing)
        #print('missing', missing)
        #print('2')
    else:
        changed = True
        #print(line)
        #print(missing)
        #print('oldac', alleleCounts)
        #print('new ac', newAC)
        #print('old alts', alts)
        #print('new alts', newAlts)
        #print('new alts ref', newAltsRef)
        #print('only missing', onlyMissing)
        #print('missing', missing)
        if missing == True:
            newAlts.append('N')
        line[4] = ','.join(newAlts)
        for a in range(len(newAC)):
            newAC[a] = str(int(int(newAC[a])/2))
        line[7] = 'AC='+ ','.join(newAC)+';'+line[7][1]
        
        #print(line)
        #print(newAC)
        #print(3)
    #print('\t'.join(line))
    #test = ['A','C']
    #print('test',','.join(test))

    #print('only', only)
    #print('both', both)
    #print(line)
    #print(line[1])
    #print('missingOnly', onlyMissing)
    #print('missing', missing)
    #print('variantOnly', variantOnly)
    return newAlts, newAltsRef, newAC, line,variantsPresent, onlyMissing, changed

def processIndels(line,newAlts, newAltsRef):
    print(line[1])
    ref = line[3]
    for a in newAlts:
        if len(ref) == len(a):
            print('snp?')
            snps = []
            for i in range(len(ref)):
                if ref[i] != a[i]:
                    snps.append(i)
            print('snps',snps)

        elif len(ref) > len(a):
            print('del', ref, a)
        else:
            print('ins', ref, a)

    '''
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
    '''

def processSNPs(line, newAlts, newAltsRef):
    missingperline = 0
    for i in range(9,len(line)):
        
        gt = line[i].split('/')
        #when written as diploid the gt should have matching haplotype
        assert gt[0] == gt[1]
        
        #if a sample has missing data, determine if 'N' needs to be added to alt alleles and change genotype to the index for 'N'
        if line[i] == './.':
            missingperline += 1
            #check this!!!!
            #if ref is 0 then all alts are 1 based encoding?
            #line[i] = str(line[4].split(',').index('N')+1)
            line[i] = str(newAlts.index('N')+1)
            
        #if no data is missing for sample, change genotype to haploid   
        else:
            if gt[0] == '0':
                line[i] = gt[0]
            else:
                oldInd = gt[0]
                #print(oldInd, type(oldInd))
                #print(newAltsRef)
                newInd = str(newAltsRef.index(oldInd)+1)
                line[i] = newInd
                #print('old',oldInd)
                #print('newInd', newInd)
    return line, missingperline

#open input for reading and open output for writing 
with gz.open(inFile, 'rb') as f:
    with open(outFile, 'w') as of:
        toProcess = 0
        totalLines = 0
        deleted = 0
        indel = 0
        snp = 0
        #iterate through input
        for line in f:
            
            #convert from bytes to str 
            line = str(line, 'UTF-8')
            
            #skip through header, but write to new vcf
            if line.startswith('#'):
                of.write(line)
                continue
            totalLines += 1
            #print(line)
            #for each line, use processAlleles to determine if the line has any relevant variants or missing data
            data = processAlleles(line)
            if data != None:
                newAlts, newAltsRef, newAC, line,variantsPresent, onlyMissing, changed = data
                #print(newAlts, newAltsRef, newAC, line, variantsPresent, onlyMissing)
                toProcess += 1
            else:
                deleted += 1
                continue

            posIndel = False 
            for a in newAlts:
                if len(a) > 1:
                    posIndel = True
                    break
            
            if posIndel == True or len(line[3]) > 1:
                print('indel? ')
                processIndels(line, newAlts, newAltsRef)
                indel += 1
            
            else:
                snp += 1
                line, missingperline = processSNPs(line, newAlts, newAltsRef)
                #print(line, missingperline)

            #note that AC does not count the number of N alleles per line 
            #print('missingperline', missingperline)
            #print('\t'.join(line))
            of.write('\t'.join(line)+'\n')

        print('processed', toProcess)
        print('deleted', deleted)
        print('indel', indel)
        print('snp', snp)
        print('total', totalLines)



#old code
'''
    missing = False
    missingperline = 0
    possibleIndel = False
    for a in alts:
        if len(a) > 1:
            possibleIndel = True
            break

    if possibleIndel == True or len(line[3]) > 1:
        print('indel?')
        #processIndels(line, alleleCounts, alts, variantPresent, onlyMissing, of)
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
'''

#this is the process alleles functions
'''
line = line.strip().split('\t')
#note that alleleCounts are 2X the amount bc it assumes diploid
alleleCounts = line[7].split(';')[0][3:].split(',')
#print('AC',alleleCounts)
alts = line[4].split(',')
#print('alts',alts)
#make sure there is a count for every allele
assert len(alts) == len(alleleCounts)


#iterate through alleleCounts to determine if there are variants for any alleles at the position

variantsPresent = False

newAlts = []
newAltsRef = []
newAC = []
for i in range(len(alleleCounts)):
    #print(alleleCounts[i])
    if int(alleleCounts[i]) > 0:
        #make sure ref is 1 indexed
        newAltsRef.append(i+1)
        newAlts.append(alts[i])
        newAC.append(alleleCounts[i])
        variantsPresent = True
#print('new ac', newAC, 'newAltsRef', newAltsRef, 'newAlts', newAlts)
#print('variants?',variantsPresent)
assert len(newAltsRef) == len(newAlts)
assert len(newAlts) == len(newAC)

onlyMissing = False
# if there are no variants, check for missing data before deleting 
if variantsPresent==False:
    if './.' not in line:
        #there is no missing data or variants, this line can be ignored
        deleted += 1
        #print('deleted')
        continue
    else:
        #there is missing data in this line and it must be processed
        missingData += 1
        onlyMissing = True
else:
    detected += 1
missing = False
missingperline = 0
toProcess += 1
print('lines to process', line, toProcess)

print('missingOnly', onlyMissing)
#line[7].split(';')[0] = 'AC='+ ','.join(newAC)

#print('new',line)
'''
#this is the process snps function
'''
missingperline = 0
for i in range(9,len(line)):
    
    gt = line[i].split('/')
    #when written as diploid the gt should have matching haplotype
    assert gt[0] == gt[1]
    
    #if a sample has missing data, determine if 'N' needs to be added to alt alleles and change genotype to the index for 'N'
    if line[i] == './.':
        missingperline += 1
        #check this!!!!
        #if ref is 0 then all alts are 1 based encoding?
        #line[i] = str(line[4].split(',').index('N')+1)
        line[i] = str(newAlts.index('N')+1)
        
    #if no data is missing for sample, change genotype to haploid   
    else:
        if gt[0] == '0':
            line[i] = gt[0]
        else:
            oldInd = gt[0]
            #print(oldInd, type(oldInd))
            #print(newAltsRef)
            newInd = str(newAltsRef.index(oldInd)+1)
            line[i] = newInd
            #print('old',oldInd)
            #print('newInd', newInd)
            '''