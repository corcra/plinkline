# The purpose of this script is to "translate" ped files. We take a phased ped file (so the first of every pair of bases in a SNP corresponds to the same chromosome, etc.), separate the two chromosomes, discard unwanted SNPs (where fox is not homozygous), translate the rest of the bases into 0 for derived and 1 for ancestral, as dictated by the "ihsmap" file (basically a modified fox .map file)

import sys
import os
import glob

def make_transfile(chro,ihsmap_folder_path):                        #this is a dictionary whose keys are 'line locations in a ped file' and values are 'derived and ancestral alleles at those sites'
    ihsmap_file_path = ihsmap_folder_path+'chr'+str(chro)+'.ihsmap'
    trans_dict = dict()
    for line in open(ihsmap_file_path,'r'):
        derived = line.split()[3]
        ancestral = line.split()[4]
        location = int(line.split()[5])
        trans_dict[location] = (derived,ancestral)
    return trans_dict

def trans_line(line,trans_dict,chro):                               # translates single line
    sex = int(line.split()[4])                                      # recall, males only have one chrX!
    translated1 = []
    translated2 = []
    linelen = len(line.split())
    for loc in range(6,linelen,2):
        base = line.split()[loc]
        nextbase = line.split()[loc+1]
        transloc = (loc-6)/2                                        # this line depends on the format of the map file
        if transloc in trans_dict:
            derived = trans_dict[transloc][0]
            ancestral = trans_dict[transloc][1]
            if base == ancestral:
                newbase = '1'
            else:
                newbase = '0'
            if nextbase == ancestral:
                newnextbase = '1'
            else:
                newnextbase = '0'
# here is a neater way of doing the above:
# newbase = (base==derived)*'0'+(base==ancestral)*'1'
# newnextbase = (nextbase==derived)*'0'+(nextbase==ancestral)*'1'
            if len(newbase)!=1 or len(newnextbase)!=1:
                print loc,'Got a problem here: apparently these bases are neither derived nor ancestral?', base, nextbase, '(derived/anc:',derived,ancestral+')'
                sys.exit('Intolerable!')
            translated1.append(newbase)
            translated2.append(newnextbase)                         # technically speaking, we don't need to do this for males on chrX, since we throw this away anyway
    if sex == 1 and chro==39:                                       # males only contribute one haplotype on chrX!
        out = ' '.join(translated1)+'\n'
    else:
        out = ' '.join(translated1)+'\n'+' '.join(translated2)+'\n'
    return out

def translate(arguments):
    chr_list = arguments[0]
    pop_name = arguments[1]
    ped_folder_path = arguments[2]
    ihshap_folder_path = arguments[3]
    ihsmap_folder_path = arguments[4]

    for chro in chr_list:
        pop_ped_path = ped_folder_path+'chr'+str(chro)+'/'+pop_name+'.chr'+str(chro)+'.ped'
        pop_ihshap_path = ihshap_folder_path+'chr'+str(chro)+'/'+pop_name+'.chr'+str(chro)+'.ihshap'

        if not os.path.exists(pop_ped_path):
            print 'No ped file exists for '+pop_name+' (chr'+str(chro)+') - skipping!'
        elif os.path.exists(pop_ihshap_path):
            print 'ihshap already exists for '+pop_name+' (chr'+str(chro)+') - skipping!'
        else:
            print 'Translating chr'+str(chro)+' of '+pop_name
            trans_dict = make_transfile(chro,ihsmap_folder_path)
            out_file = open(pop_ihshap_path,'w')
            ped_file = open(pop_ped_path,'r')
#i=0
            for line in ped_file:
#               i+=1
#               if i>30:                                            #this is how we make it first30
#                   break
#               else:
                translated = trans_line(line,trans_dict,chro)
                out_file.write(translated)
            out_file.close()

def run_translation(chr_lists,pop_names,ped_folder_path,ihshap_folder_path,ihsmap_folder_path,p,cores):
    for pop_name in pop_names:
        if len(glob.glob(ihshap_folder_path+'chr1/'+pop_name+'*'))==0:      # the assumption is that if chr1 has been translated, so too have all of the rest (since chr1 usually goes last)
            print 'Translating population',pop_name
            arg_list = zip(chr_lists,[pop_name]*cores,[ped_folder_path]*cores,[ihshap_folder_path]*cores,[ihsmap_folder_path]*cores)
            p.map(translate,arg_list)
        else:
            print 'Have already translated',pop_name                # this is very clearly not robust, but it should work well enough...
