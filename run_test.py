# This basically just calls the Pritchard scripts!

import glob
import os
from subprocess import call

def run_ihs(arguments):
    chr_list = arguments[0]
    pop_name = arguments[1]
    statistics_folder_path = arguments[2]
    ihshap_folder_path = arguments[3]
    ihsmap_folder_path = arguments[4]
    ihs_path = arguments[5]
    test_type = arguments[6]

    already_ihs = []

    for chro in chr_list:
        out_file_path = statistics_folder_path+test_type+'/chr'+str(chro)+'/'+pop_name+'.chr'+str(chro)
        if len(glob.glob(out_file_path+'*'))!=0:            # if the data already exists for this chromosome, don't do it again
            already_ihs.append(str(chro))
        else:
            print 'Running ihs on chr'+str(chro)
            pop_path = ihshap_folder_path+'chr'+str(chro)+'/'+pop_name+'.chr'+str(chro)+'.ihshap'
            ihsmap_path = ihsmap_folder_path+'chr'+str(chro)+'.ihsmap'
            if os.path.exists(pop_path):
                if os.path.exists(ihsmap_path):
                    script = ihs_path+' '+ihsmap_path+' '+pop_path+' > '+out_file_path
                    call(script,shell=True)
                else:
                    print 'Missing file:', ihsmap_path
            else:
                print 'Missing file:', pop_path

    if len(already_ihs)>0:
        print 'Already had iHS for', pop_name, 'on chrs '+', '.join(already_ihs)

def run_xpehh(arguments):
    chr_list = arguments[0]
    pop_1_name = arguments[1]
    pop_2_name = arguments[2]
    statistics_folder_path = arguments[3]
    ihshap_folder_path = arguments[4]
    ihsmap_folder_path = arguments[5]
    xpehh_path = arguments[6]
    test_type = arguments[7]

    extract_xpehh = str()
    already_xpehh = []

    for chro in chr_list:
        out_file_path = statistics_folder_path+test_type+'/chr'+str(chro)+'/'+pop_1_name+'-'+pop_2_name+'.chr'+str(chro)
        if len(glob.glob(out_file_path+'*'))!=0:            # if the data already exists for this chromosome, don't do it again
            already_xpehh.append(str(chro))
        else:
            print 'Running XP-EHH on chr'+str(chro)
            pop_1_path = ihshap_folder_path+'chr'+str(chro)+'/'+pop_1_name+'.chr'+str(chro)+'.ihshap'
            pop_2_path = ihshap_folder_path+'chr'+str(chro)+'/'+pop_2_name+'.chr'+str(chro)+'.ihshap'
            ihsmap_path = ihsmap_folder_path+'chr'+str(chro)+'.ihsmap'

            if os.path.exists(pop_1_path):
                if os.path.exists(pop_2_path):
                    if os.path.exists(ihsmap_path):
                        script = xpehh_path+' -m '+ihsmap_path+' -h '+pop_1_path+' '+pop_2_path+' > '+out_file_path
                        call(script,shell=True)
                    else:
                        print 'Missing file:', ihsmap_path
                else:
                    print 'Missing file:', pop_2_path
            else:
                print 'Missing file:', pop_1_path

    if len(already_xpehh)>0:
        print 'Already had XP-EHH for', pop_1_name,'and', pop_2_name, 'on chrs '+', '.join(already_xpehh)

def run_the_test(pop_names,chr_lists,test_type,statistics_folder_path,ihshap_folder_path,ihsmap_folder_path,ihs_path,xpehh_path,p,cores):
    if test_type == 'ihs':
        pop_name = pop_names[0]
        print 'Running iHS on '+pop_name+'!'
        arg_list = zip(chr_lists,[pop_name]*cores,[statistics_folder_path]*cores,[ihshap_folder_path]*cores,[ihsmap_folder_path]*cores,[ihs_path]*cores,[test_type]*cores)                              # the reason this is a horrible mess is that p.map only allows you to pass a single argument, so I zip them together into a tuple of sorts here
        p.map(run_ihs,arg_list)
    else:
        pop_1_name = pop_names[0]
        pop_2_name = pop_names[1]
        print 'Running XP-EHH on '+pop_1_name+' and '+pop_2_name+'!'
        arg_list = zip(chr_lists,[pop_1_name]*cores,[pop_2_name]*cores,[statistics_folder_path]*cores,[ihshap_folder_path]*cores,[ihsmap_folder_path]*cores,[xpehh_path]*cores,[test_type]*cores)
        p.map(run_xpehh,arg_list)
