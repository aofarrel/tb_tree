import os
import argparse
import gzip
import sys

#parser = argparse.ArgumentParser()
#parser.add_argument('-v', '--megaVCF', required=True, type=str,help='path to giant VCF to be chunked')
#parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
#parser.add_argument('-c', '--chunkSize', type=int, help='number of samples included in file', default=1000)
#parser.add_argument('-t', '--threads', type=int, help='number of threads (where applicable)', default=1)
#parser.add_argument('-sd', '--scriptsDirectory', type=str, required=True, help='path to directory where scripts are')


'''
THINGS TO NOTICE:
Unix `cut` command 1-indexes file columns

I will be using linux's column ordering system to track columns that are processed 

if one wants to see if a sample ('column') has been processed, it will be the same in the processed file name and the orignial linux ordering

file names are inclusive of all numbers 

diff file lists 1st position and the full lenght of the string w this char, this may need to change becasue the lenght value is inclusive of the first sight 
'''

#args = parser.parse_args()
#vcf = args.megaVCF
#chunk = args.chunkSize
#wd = args.workingDirectory

file = sys.argv[1]

#makes sure input path wont cause error
'''
if wd[-1] != '/':
    print('add a final slash')
    wd = wd+'/'
sd = args.scriptsDirectory
if sd[-1] != '/':
    print('add a final slash')
    sd = sd+'/'

t = args.threads
'''

def squish(file):
    with open(file) as f:
        with open(f'{file}squish', 'w') as o:

            prev = None
            for line in f:
                if not line.startswith('>'):

                    #process line
                    line = line.strip().split()
                    print(line, prev)
                    
                    if prev != None:
                        print('line', line)
                        print('prev', prev)
                        if prev[0] == line[0]:
                            print('prev', prev, 'line', line)
                            if int(prev[1]) == int(line[1])-int(prev[2]):
                                print('prev', prev, 'now', line)
                                prev[2] = str(int(prev[2])+int(line[2]))
                                print('new prev', prev)
                            else:
                                o.write('\t'.join(prev)+'\n') 
                                prev = line
                        else:
                            o.write('\t'.join(prev)+'\n') 
                            prev = line
                        
                    
                    #the first line of file becomes prev variable 
                    else:
                        prev = line
                        print('first prev!', prev)
                        
                        
            

                #write header to new file
                else:
                    o.write(line)
                    



squish(file)