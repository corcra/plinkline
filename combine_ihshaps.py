#!/usr/bin/python

#this will take a list of breeds TO COMBINE and combine their ihshap files

import sys
from subprocess import call
import locations

# argv[1]: file with populations (eg corresponding to existing ihshap files) to combine
# argv[2]: number of haplotypes to take, per population (for allosomes this equates to half the number of individuals, of course)

if len(sys.argv)<2:
    sys.exit('Usage: combine_ihshaps.py FILENAME num_hap\nIf num_hap = 0, all individuals/haplotypes in the subpopulations are taken, else num_hap from each.')

chr_list = range(1,40)

ihshap_path = locations.ihshap_folder_path
out_path = locations.ihshap_folder_path
joblist_name = sys.argv[1]

per_pop = sys.argv[2]

if per_pop == '0':
    do_line = 'cat'
else:
    do_line = 'head -'+per_pop

pops_to_combine = [line.split()[0] for line in open(joblist_name)]

for chro in chr_list:
    chr_path = ihshap_path+'chr'+str(chro)+'/'
    pops_file_paths = [chr_path+pop+'.chr'+str(chro)+'.ihshap' for pop in pops_to_combine]
    pops_joined = ' '.join(pops_file_paths)
    string = 'for f in '+pops_joined+'; do '+do_line+' $f >> '+out_path+'chr'+str(chro)+'/'+joblist_name+'.chr'+str(chro)+'.ihshap; done'
    call(string,shell=True)
