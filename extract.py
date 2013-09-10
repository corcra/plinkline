# This script will be called by master_plinkline.py
# Purpose: given a large ped file with many individuals, extract ones we want (known by their ID) and put them into a new ped file (having the same name as the name of the file containing the IDs).

import os
import sys
import glob

def get_ids(pop_names,keeplists_folder_path):
    id_pop = dict()                                                         # key is ID, value is population name
    for pop_name in pop_names:
        pop_keeplist_file_path = keeplists_folder_path+pop_name+'.idlist'
        if not os.path.exists(pop_keeplist_file_path):
            print 'No ID list for', pop_name, '- cannot extract!'
        else:
            keeplist_file = open(pop_keeplist_file_path,'r')
            for line in keeplist_file:
                ID = line.split()[0]
                id_pop[ID] = pop_name
            keeplist_file.close()
    return id_pop

def split_ped(arguments):
    chr_list = arguments[0]
    id_pop = arguments[1]
    pop_files = arguments[2]
    ped_folder_path = arguments[3]
    master_ped_name = arguments[4]

    for chro in chr_list:
        ped_path = ped_folder_path+master_ped_name+'.chr'+str(chro)+'.ped'
        if os.path.exists(ped_path):
            print 'Extracting chr'+str(chro)+'...'
            out_ped_files = dict()
            pop_counts = dict()                                             # for the purpose of recalling how many of each population we recorded
            for pop_name in set(id_pop.values()):
                pop_ped_path = ped_folder_path+'chr'+str(chro)+'/'+pop_name+'.chr'+str(chro)+'.ped'
                out_ped_files[pop_name] = open(pop_ped_path,'w')            # dictionary of file objects to write to
                pop_counts[pop_name] = 0

            for ped_line in open(ped_path,'r'):
                ID = ped_line.split()[0]
                if ID in id_pop:
                    pop_name = id_pop[ID]
                    out_ped_files[pop_name].write(ped_line)
                    pop_counts[pop_name]+=1

            for pop_name in set(id_pop.values()):
                print '\tExtracted '+str(pop_counts[pop_name])+' '+pop_name+' individuals on chr '+str(chro)+'.'
                out_ped_files[pop_name].close()

        else:
            print 'Can\'t find master ped file on chr'+str(chro)+' - cannot extract!'
            print ped_path



# --- the main 'loop' is here ! --- #
def run_extraction(pop_names,keeplists_folder_path,chr_lists,ped_folder_path,master_ped_name,p,cores):
    to_extract = []
    for pop_name in pop_names:                                              # pop_files is a dict of pop_name:pop_file_stem
        if len(glob.glob(ped_folder_path+'chr1/'+pop_name+'*'))==0:		    # chr1 is actually one of the last chromosomes that will have been extracted, so we check to see if it exists. If it doesn't we assume none do, which is a strong assumption!
            if len(glob.glob(keeplists_folder_path+pop_name+'.idlist'))!=0:
                to_extract.append(pop_name)
            else:
                print 'Have no ped file for', pop_name, 'but missing ID file - cannot extract!'
                print ped_folder_path+'chr1/'+pop_name                      # probably something has gone wrong if you're seeing this...
                print keeplists_folder_path+pop_name+'.idlist'

    if len(to_extract)==0:
        print 'No populations to extract!'
    else:
        print 'Creating individual ped files for '+', '.join(to_extract)+'...'
        id_pop = get_ids(to_extract,keeplists_folder_path)                  # id_pop is a dict of id:pop_name (as decided by the name of the file the ID came from)
        if len(id_pop)!=0:
            arg_list = zip(chr_lists,[id_pop]*cores,[pop_names]*cores,[ped_folder_path]*cores,[master_ped_name]*cores)
            p.map(split_ped,arg_list)
