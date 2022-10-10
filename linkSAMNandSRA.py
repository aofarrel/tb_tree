import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-sd', '--sradir', required=True)
#parser.add_argument('-d', '--workingDirectory', required=True, type=str, help='directory for all outputs (make sure this directory will have enough space!!!!)')
parser.add_argument('-o', '--conversion_file', required=True)


args = parser.parse_args()
sd = args.sradir
output = args.conversion_file

#makes sure input path wont cause error
if sd[-1] != '/':
    print('add a final slash')
    sd = sd+'/'

samples = {}
for file in os.listdir(sd):
    #print(file)
    if file.endswith('json'):
        SAMN = file[:-5]
        if SAMN not in samples:
            samples[SAMN] = {'samn-last-update': None, 'SRAs': None, 'geo-location': None}
        #else:
        #    print('json',samples[SAMN])
        #geo = sys.stdin(os.system(f'cat {sd}{file} | grep geo'))
        #print('geo', geo)
        
        with open(f'{sd}{file}') as f:
            for line in f:
                #print(line)
                if 'geo_loc_name' in line and samples[SAMN]['geo-location'] == None:
                    #print(line)
                    line = line.strip().split('": ')
                    #print(line[1][1:-2])
                    #samples[SAMN]['geo-location'] == 
                    samples[SAMN]['geo-location'] = line[1][1:-2]

                if 'ENA-LAST-UPDATE' in line and samples[SAMN]['samn-last-update'] == None:
                    line = line.strip().split('": ')
                    #print(SAMN, line[1][1:-1])
                    samples[SAMN]['samn-last-update'] =  line[1][1:-1]
                    break
        #print(SAMN, samples[SAMN])
            
    elif file.endswith('accession'):
        #print(file)
        SAMN = file[:-9]
        if SAMN not in samples:
            samples[SAMN] = {'samn-last-update': None, 'SRAs': None, 'geo-location': None}
        #else:
        #    print('accession', samples[SAMN])

        with open(f'{sd}{file}') as f:
            for line in f:
                print(file, line)
                if samples[SAMN]['SRAs'] == None:
                    print(line.strip())
                    samples[SAMN]['SRAs'] = [line.strip()]
                else:
                    samples[SAMN]['SRAs'].append(line.strip())

with open(output,'w') as o:
    o.write('SAMN_ID\tsamn_last_update\tSRAs\tgeolocation\n')
    for l in samples:
        for key in samples[l]:
            line = [samples[l][key] for key in samples[l]]
        o.write(f'{l}\t{line[0]}\t{",".join(line[1])}\t{line[2]}\n')
            
        #print(l,*(samples[l][key] for key in samples[l]), sep='\t')

        


                #if 'geo_loc_name' in line and samples[SAMN]['geo-location'] != None:
                #    print('already has geo?',line)
                    

                    

                #print(file, line)


