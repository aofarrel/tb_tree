import argparse
import gzip as gz
import copy

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcfFile', required=True)
#parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-o', '--outFile', required=True)
args = parser.parse_args()
inFile=args.vcfFile
outFile=args.outFile

#to do:
# change missing 'N' to matching len of reference 
# edge cases 

def processLines(line):
    line = line.strip().split('\t')
    alleleCounts = line[7].split(';')[0][3:].split(',')
    alts = line[4].split(',')
    ref = line[3]

    #make sure there is a count for every allele
    #there should be no alleles or counts for missing data
    assert len(alts) == len(alleleCounts)
    
    #track if there are even variants (this was already done in processVCF so shouldnt catch anything)
    variantsPresent = False
    newAlts = []
    newAltsRef = []
    newAC = []

    #already asserted that len(alts) is the same so this value should be consistent 
    for i in range(len(alleleCounts)):
        #print(alleleCounts)

        #essentially the same thing as in processVCF
        if int(alleleCounts[i]) > 0:
            #make sure alts are 1 indexed
            newAltsRef.append(str(i+1))
            newAlts.append(alts[i])
            newAC.append(alleleCounts[i])
            variantsPresent = True

    #make sure all new stats are equivalent 
    assert len(newAltsRef) == len(newAlts)
    assert len(newAlts) == len(newAC)

    #print('newAlts', newAlts)
    #print('newAltsRef', newAltsRef)

    onlyMissing = False
    missing = False

    #probably dont need this 
    variantOnly = False

    # if there are no variants, check for missing data before deleting 
    #this is redundant but probably too complicated to easily delete
    if variantsPresent==False:
        if './.' not in line:
            #there is no missing data or variants, this line can be ignored
            return None
        else:
            #there is missing data in this line and it must be processed
            onlyMissing = True
    
    #variants present
    else:
        #determine if there are variants and missing data
        if './.' in line:
            missing = True
        else:
            variantOnly = True
    
    #there should either be missing data or variants (allows for both as well)
    assert onlyMissing == True or len(newAC) > 0
    
    #change allele counts (counts changed from dip to hap)
    #if greater efficiency is needed, dip->hap is not strictly necessary
    line[7] = line[7].split(';')
    #print('line[7]',line[7])
    changed = False

    #what are these?
    one = 0
    two = 0
    three = 0

    #the first conditional accounts for all lines where onlyMissing = True
    if len(newAC) == 0:
        #print(line[1])

        #add string N to represent missing data
        #track this later to see if it confuses usher
        newAlts.append('N'*len(ref))
        #actually update the line
        line[4] = 'N'*len(ref)

        #do not create an allele count for missing data
        #add an back in 
        line[7] = 'AC='+';'+line[7][1]
        #track if the line was updated
        changed = True
        #print('1')

    #this conditional note's that no alt alleles were removed from the position and updates format wo changing anything else
    #
    elif alleleCounts == newAC:
        #determine if there is missing data and add to alts and line
        #print('no change, update format') 
        if missing == True:
            newAlts.append('N'*len(ref))
            line[4] = ','.join(newAlts)
            #print('there is missing in this!!!', line)
        
        #is this not allele counts? fix
        #print('allele counts', alleleCounts)
        ac = line[7][0][3:].split(',')
        #print('ac', ac)

        #probably dont need this anymore? i control the ac now
        ####THIS IS WEIRD FIGURE OUT WHAT IT"S DOING 
        #for a in range(len(ac)):
        #    ac[a] = str(int(int(ac[a])/2))
        #print('ac after', ac)

        line[7] = 'AC='+ ','.join(ac)+';'+line[7][1]
        #print(line[7])
    
    #this conditional is if there are alleles but less than the og file 
    else:
        #print('new allele list ', line)
        #do i need to track changed? just means changed from previous file
        changed = True
        
        #add alt for missing, make sure this is ok w usher later
        if missing == True:
            newAlts.append('N'*len(ref))
        
        #####come back here to check 
        #print('new alts', newAlts)
        #print('newAC', newAC)

        line[4] = ','.join(newAlts)
        #line[7]
        #print(line)
        #for a in range(len(newAC)):
        #    print('a', a)
        #    print(int(int(newAC[a])/2))
        #    newAC[a] = str(int(int(newAC[a])/2))
        #print('new ac?',newAC)
        line[7] = 'AC='+ ','.join(newAC)+';'+line[7][1]
        #print('updated?', line)
        
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


def processIndels(line,newAlts, newAltsRef, newAC):
    addLines = []
    removeAlt = []
    print('for indels', line)
    print('newAlts', newAlts)
    print('newAltsRef', newAltsRef)
    print('newAC', newAC)
    #snp = False
    #print('work?')
    ref = line[3]
    #print('position',line[1])
    #print('preprocess',line)
    #print(newAlts)
    #print(newAltsRef)
    #newLine = []
    line, missingperline = processAlleles(line, newAlts, newAltsRef)
    dontAdd = False
    #print('postprocess',line)
    
    for a in range(len(newAlts)):

        if len(ref) == 1:
            print(ref,newAlts, newAlts[a])
            if len(ref) == len(newAlts[a]):
                assert newAlts[a] == 'N'
            else:
                assert len(ref) < len(newAlts[a])
            #the simplest case: insertion from one base pair to many and no other alleles... simply delete line from file
            #misrepresents these sites as reference but that is currently how usher does it 
            #if there are no other alts in the site there is no reason to store it even if this misrepresents it
            if len(newAlts) == 1:
                #print('delete', line)
                return None
            
            #if there is more than one alt, i don't think i should change gt back to reference for insertions bc the variant is 
            #still in the file), replace the alt with a line of Ns so we know that it's an indel but also that it's been processed
            elif len(newAlts) > 1:
            
                #replacing the alt with a string of N's represents that it's an indel wo removing it from vcf and misrepresenting
                #site as reference genotype


                newAlts[a] = 'N'*len(newAlts[a])
                #line[4] = ','.join(newAlts)
                #change alts to updated newAlts
                changeAlt = line[4].split(',')
                changeAlt[a] = 'N'*len(newAlts[a])
                line[4] = ','.join(changeAlt)
                #note that nothing here is changed except how the indel is represented (we are changing alt so we can see that its 
                #been processed) usher will see this site as the reference bc it ignores indels but we will know there is a variant here

        #the less simple cases 
        else:
            #if len(ref) is less than allele but greater than 1, there is likely something strange about this site, better to mask it than try to 
            #process
            #MAKE SURE USHER CAN HANDLE MY FORMAT
            if len(ref) < len(newAlts[a]):
                #if there is only one alt, and we are masking the alt, just delete the line entirely (it's not quite accurate as a VCF but Usher 
                #wont know the diff)
                if len(newAlts) == 1:
                    #print('delete', line)
                    return None
                #if there is more than one alt, we are going to mask this insertion in the line of alts, bc we are keeping the site and showing that it
                #has been prpcessed (this approach may need to be modified eventually)
                else:
                    newAlts[a] = 'N'*len(newAlts[a])
                    changeAlt = line[4].split(',')
                    changeAlt[a] = 'N'*len(newAlts[a])
                    line[4] = ','.join(changeAlt)
                    #line[4] = ','.join(newAlts)

            #deletions need subsequent positions masked as well
            #since subsequent positions will be added, need a new line for everysingle one 
            elif len(ref) > len(newAlts[a]):
                #first modify the alt to a string of Ns
                print(a)
                newAlts[a] = 'N'*len(newAlts[a])
                changeAlt = line[4].split(',')
                changeAlt[a] = 'N'*len(newAlts[a])
                print(changeAlt)
                line[4] = ','.join(changeAlt)
                #line[4] = ','.join(newAlts)
                #don't need the intial deletion if it is the only variant in the position
                #misrepresents indel as reference but usher cant tell 
                #THIS MAY NEED TO BE MODIFIED, SOME DELETIONS MAY NOT LINE UP FULLY KEEP AN EYE OUT
                if len(newAlts) == 1:
                    dontAdd = True
                
                #create new lines for all deleted sites from end of alt to end of ref
                end = len(newAlts[a])
                extra = ref[end:]
                #the line has been modified, we want to make a copy so that we are 
                #not accidentally creating and updating a reference 
                newLine = copy.copy(line)
                #modify newline for each subsequent newline (they all have the same alt and samples with the genotype)
                newLine[4] = 'N'
                newLine[7] = 'AC='+';'+line[7][1]
                for j in range(9, len(newLine)):
                    #print(j, newLine[j], a)
                    if newLine[j] == str(a+1):
                        #print('inline',j, a, a+1)
                        newLine[j] = '1'
                    else:
                        newLine[j] = '0'
                #iterate through all deleted positions and create a line for each one 
                #add lines to list of positions to be added to vcf
                #might need to merge some of these lines later
                #worth creating a generator from this or no?
                for i in range(len(extra)):
                    nLine = copy.copy(newLine)
                    nLine[1] = str(i + int(line[1])+end)
                    nLine[3] = extra[i]
                    addLines.append(nLine)
      
            #this is not an indel it is a series of close together snps 
            #need to identify all snps and add a new line for each one 
            #need to remove snp from new alts or remove line entirely if its the only alt
            elif len(ref)==len(newAlts[a]) and newAlts[a] != 'N'*len(ref):
                end = len(ref)
                snps = []
                for i in range(end):
                    if ref[i] != newAlts[a][i]:
                        snps.append(i)
                    
                #print('snps', snps)
                #create a new line that can be modified wo messing up line
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
                    #print('remove snp from this', line)
                    #print(newAlts)
                    #newAlts = newAlts.pop(a)
                    #print('updated new alts', newAlts)
                    #print('old line 4', line[4])

                    #update alts
                    rmSNP = line[4].split(',')
                    #print(rmSNP)
                    del rmSNP[a]
                    line[4] = ','.join(rmSNP)
                    #print('new line 4', line[4])

                    #update counts
                    #print('old line 7', line[7])
                    changeCT = line[7].split(';')[0][3:].split(',')
                    #print('change ct', changeCT)
                    del changeCT[a]
                    #this is hyperspecific and not the most important thing if it stops working dont mess w it too much
                    line[7] = line[7][0:3]+','.join(changeCT) + ';' + line[7].split(';')[1]
                    
                    #print('after',changeCT)
                    #print('new line 7', line[7])
                    #print('old', line)
                    #remove alt from gts as well
                    for i in range(9, len(line)):
                        if line[i] == str(a+1):
                            line[i] = '0'
                    #print('new',line)
                    '''

    if dontAdd == False:
        if len(removeAlt) > 0:
            for a in removeAlt:
                #print('remove snp from this', line)
                #print(newAlts)
                #newAlts = newAlts.pop(a)
                #print('updated new alts', newAlts)
                #print('old line 4', line[4])

                #update alts
                rmSNP = line[4].split(',')
                print('line[4]',rmSNP)
                print('allele', a)
                del rmSNP[a]
                line[4] = ','.join(rmSNP)
                #print('new line 4', line[4])

                #update counts
                #print('old line 7', line[7])
                changeCT = line[7].split(';')[0][3:].split(',')
                #print('change ct', changeCT)
                del changeCT[a]
                #this is hyperspecific and not the most important thing if it stops working dont mess w it too much
                line[7] = line[7][0:3]+','.join(changeCT) + ';' + line[7].split(';')[1]
                
                #print('after',changeCT)
                #print('new line 7', line[7])
                #print('old', line)
                #remove alt from gts as well
                for i in range(9, len(line)):
                    if line[i] == str(a+1):
                        line[i] = '0'
                #print('new',line)

        addLines.insert(0,line)
    lines = checkNewLines(addLines)
    return lines
    #print('line',line)
    #print('newlines', addLines)

def checkNewLines(addLines):
    positions = {}
    lines = []
    for l in range(len(addLines)):
        if int(addLines[l][1]) not in positions:
            positions[int(addLines[l][1])] = [l]
        else:
            #positions[int(addLines[l][1])][0] += 1
            positions[int(addLines[l][1])].append(l)
    for p in positions:
        print(p, positions[p])
        if len(positions[p]) > 1:
            newLine = copy.copy(addLines[positions[p][0]])
            #print('old line',addLines[positions[p][0]])
            #print(newLine)
            for l in range(1,len(positions[p])):
                #newLine
                #print('line',addLines[positions[p][l]])
                #need to make sure this can work for any lenght of alts
            
                assert newLine[1] == addLines[positions[p][l]][1]
                assert newLine[3] == addLines[positions[p][l]][3]
                alts = newLine[4].split(',')
                alAlts = addLines[positions[p][l]][4].split(',')
                print('alts',alts)
                print('alakts',alAlts)
                for i in range(len(alAlts)):
                    allele = i+1 
                    alts.append(alAlts[i])
                    #print('alts', alts)
                    for j in range(9,len(addLines[positions[p][l]])):
                        #print(j, addLines[positions[p][l]][j])
                        if addLines[positions[p][l]][j] == str(allele):
                            assert newLine[j] == '0'
                            newLine[j] = str(alts.index(alAlts[i])+1)
            newLine[4] = ','.join(alts)
            print('newline', newLine)
            lines.append(newLine)
        else:
            lines.append(addLines[positions[p][0]])


    lines = sorted(lines, key=lambda x: x[1])
    #print('indel lines', lines)
    return lines

                
    #need to check to make sure generator stays in order 
    #yield line
        

        
    '''
        if len(ref) == len(newAlts[a]):
            print('snp?')
            snps = []
            print(ref)
            print(a)
            for i in range(len(ref)):
                
                if ref[i] != newAlts[a][i]:
                    snps.append(i)
                    newLine.append(['' for j in range(len(line))])
            print(newLine)
            print(snps)
            for snp in range(len(snps)):
                print(newLine[snp])
            
        
            if len(snps) > 1:
                print(newline)
            elif len(snps) == 1:
                line[1] = str(snps[0] + int(line[1]))
                line[3] = line[3][snp[0]]
                line[4] = line[4]
            else:
                #put an actual error checker here
                print('error!!!!, this isnt a variant ')
            
            

        elif len(ref) > len(newAlts[a]):
            print('del', ref, a)
        else:
            print('ins', ref, a)
        '''

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

def processAlleles(line, newAlts, newAltsRef):
    print('for snps!!!', line[1])
    #print('new alts', newAlts)
    #print('new alts ref', newAltsRef)
    missingperline = 0
    
    ref = line[3]
    #print('ref', ref)
    for i in range(9,len(line)):
        
        gt = line[i].split('/')
        #print('gt', gt)
        #when written as diploid the gt should have matching haplotype
        assert gt[0] == gt[1]
        
        #if a sample has missing data, determine if 'N' needs to be added to alt alleles and change genotype to the index for 'N'
        if line[i] == './.':
            missingperline += 1
            #check this!!!!
            #if ref is 0 then all alts are 1 based encoding?
            #line[i] = str(line[4].split(',').index('N')+1)
            line[i] = str(newAlts.index('N'*len(ref))+1)
            #print('new line', line)
            #print('there is missing in this!!!', line)
            
        #if no data is missing for sample, change genotype to haploid   
        else:
            
            if gt[0] == '0':
                line[i] = gt[0]
            else:
                oldInd = gt[0]
                #print('need to change allele' )
                #print(oldInd, type(oldInd))
                #print(newAltsRef)
                #print(newAltsRef)
                #print(line)
                newInd = str(newAltsRef.index(oldInd)+1)
                line[i] = newInd
                #print('old',oldInd)
                #print('newInd', newInd)
                #print('updated line', line)
    return line, missingperline



#MAIN (put in main funct?)

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
            
            #for each line, use processLines to determine if the line has any relevant variants or missing data
            data = processLines(line)
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
            #this needs to wirte indels to file (keep checking for edge cases)
            if posIndel == True or len(line[3]) > 1:
                pass

                #comment out to see if code works on nonindels
                #need to fix process indels 
                '''
                print('indel? ')
                #gen = 
                indelLines = processIndels(line, newAlts, newAltsRef, newAC)
                #print(indelLines)
                #for g in gen:
                #    print(g)
                indel += 1
                if indelLines != None:
                    print(indelLines)
                    #indelLines = sorted(indelLines, key=lambda x: x[1])
                    for ln in indelLines:
                        of.write('\t'.join(ln)+'\n')
                        print('order>?',ln)
                        '''
                        
            #this will write all non-indel variants into file (i think this is working well)
            else:
                snp += 1
                line, missingperline = processAlleles(line, newAlts, newAltsRef)
                #print(line, missingperline)
                #note that AC does not count the number of N alleles per line 
                #print('mqissingperline', missingperline)
                #print('\t'.join(line))
                of.write('\t'.join(line)+'\n')

        print('processed', toProcess)
        print('deleted', deleted)
        print('indel', indel)
        print('snp', snp)
        print('total', totalLines)


