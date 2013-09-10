# This just takes some output from the Pritchard scripts and produces normalised test-statistics. XP-EHH is naively normalised across the entire genome, while iHS is normalised in frequency bins of 5%.

import sys
import os
import time
import numpy as np
from subprocess import call

# Dealing with nans...
def my_sum(array):
    total = 0
    for element in array:
        if np.isfinite(element):
            total = total + element
    return float(total)

def mean_nan(array):
    denominator = my_sum(np.isfinite(array))
    numerator = my_sum(array)
    if denominator !=0:
        result = numerator/denominator
    else:
        result = np.nan
    return result

def std_nan(array):
    mean = mean_nan(array)
    devs = []
    if np.isfinite(mean):
        for element in array:
            if not np.isfinite(element):
                devs.append(np.nan)
            else:
                dev = (element-mean)**2
                devs.append(dev)

        mean_devs = mean_nan(devs)

        if np.isfinite(mean_devs):
            return np.sqrt(mean_devs)
        else:
            return np.nan
    else:
        print 'You\'re asking me to compute a standard deviation with an undefined mean? Good luck with that.\n'
        return np.nan

def make_aggregate_file(job_name,test_type,statistics_path,aggregate_file_path):
    print 'Didn\'t find aggregate file: making one!'
    if test_type == 'ihs':
        fields = '$1,$3,$6'
    else:
        fields = '$5'
    
    for chro in range(1,40):
        unnormed_path = statistics_path+test_type+'/chr'+str(chro)+'/'+job_name+'.chr'+str(chro)
        normed_path = unnormed_path+'.norm'
        if os.path.exists(unnormed_path):
            the_path = unnormed_path
            script = 'awk \'{print '+fields+' }\' '+the_path+' >> '+aggregate_file_path
            call(script,shell=True)
        else:                   # this shouldn't ever happen
            print 'Cannot find statistic file for',job_name,'on chr'+str(chro)

#this produces a dictionary with snpid: [mean,std] for normalising!
def get_norm_std(job_name,test_type,statistics_path):
    aggregate_file_path = statistics_path+test_type+'/'+job_name+'.all.'+test_type
    if not os.path.exists(aggregate_file_path):
        make_aggregate_file(job_name,test_type,statistics_path,aggregate_file_path)
        print 'Aggregate file made!'
    else:
        print 'Found aggregate file!'
    num_wins = float()
    snps_to_moments = dict()
    if test_type == 'ihs':
        window_size = 0.05
        included = 0
        master_list = []
        freq_lists = dict()
        freq_mean_stds = dict()

        for i in range(20):
            freq_lists[i] = []

        for line in open(aggregate_file_path,'r'):
            snpid = line.split()[0]
            freq = float(line.split()[1])
            value = line.split()[2]
            if value=='nan':
                print 'value is NaN!'
                print line
                value = np.nan
            else:
                value = float(value)
            master_list.append([snpid,freq,value])

        for entry in master_list:
            freq = float(entry[1])
            window_max = 0.0
            window_min = 0.0
            for i in range(20):
                window_min = window_max
                window_max = window_max + window_size
                if window_min <= freq < window_max:
                    value = entry[2]
                    freq_lists[i].append(value)
                    entry.append(i)
                    included+=1

        for i in range(1,19):
            if len(freq_lists[i])>0:
                mean = mean_nan(np.array(freq_lists[i]))
                if not np.isfinite(mean):
                    print 'mean undefined for freq window', i
                    freq_mean_stds[i] = [np.nan,np.nan]
                else:
                    std = std_nan(np.array(freq_lists[i]))
                    freq_mean_stds[i] = [mean,std]
            else:
                print 'There are no values in frequency window', i

        for entry in master_list:
            snpid = entry[0]
            window = entry[3]
            mean = freq_mean_stds[window][0]
            std = freq_mean_stds[window][1]
            if snpid in snps_to_moments:
                print 'snpid already included? what?',snpid,str(chro)
            snps_to_moments[snpid] = [mean,std]

        if len(snps_to_moments) != included:
            print len(snps_to_moments), included
            sys.exit('some problem here')
    else:
        values = [float(line) for line in open(aggregate_file_path)]
        val_array = np.array(values)
        mean = mean_nan(val_array)
        std = std_nan(val_array)
        #sort this out later
        snps_to_moments[''] = [mean,std]

    return snps_to_moments

def normalise_chr(chr_list, job_name, test_type, s_to_m, statistics_path):
    for chro in chr_list:
        file_path = statistics_path+test_type+'/chr'+str(chro)+'/'+job_name+'.chr'+str(chro)

        if os.path.exists(file_path+'.norm'):
            print 'Chromosome '+str(chro)+' is already normalised - skipping!'
        
        elif os.path.exists(file_path):
            print 'Normalising chromosome',str(chro)
            out_file_path = file_path+'.norm'
            out_file = open(out_file_path,'w')
            if test_type == 'ihs':
                for line in open(file_path):
                    snpid = line.split()[0]
                    mean = float(s_to_m[snpid][0])
                    std = float(s_to_m[snpid][1])
                    value = float(line.split()[5])
                    normed_value = (value-mean)/std
                    out_file.write(line.strip()+' '+str(normed_value)+'\n')
            else:
                mean = s_to_m[''][0]
                std = s_to_m[''][1]
                for line in open(file_path):
                    value = float(line.split()[4])
                    normed_value = (value-mean)/std
                    out_file.write(line.strip()+' '+str(normed_value)+'\n')
        # this should not be necessary, shouldn't get this far if any are missing
        else:
            print 'Cannot find', file_path,'- cannot normalise!'

def run_normalisation(test_data_name,test_type,statistics_path):
    job_name = test_data_name
    if test_type == 'xpehh':
        print 'XP-EHH: performing genome-wide normalisation without frequency bins on '+job_name
    elif test_type == 'ihs':
        print 'iHS: performing genome-wide normalisation with 5% frequency bins on '+job_name

    chr_all = range(1,40)
    already_normed = set()
    for chro in chr_all:
        file_path = statistics_path+test_type+'/chr'+str(chro)+'/'+job_name+'.chr'+str(chro)
        if os.path.exists(file_path):
            if os.path.exists(file_path+'.norm'):
                print 'chr'+str(chro),'is already normalised - skipping.'
                already_normed.add(chro)
        else:
            print 'ERROR: Missing file', file_path, '- aborting.' 
            return
    
    for chro in already_normed:
        chr_all.remove(chro)

    if len(chr_all)!=0:
        s_to_m = get_norm_std(job_name,test_type,statistics_path)
        normalise_chr(chr_all,job_name,test_type,s_to_m,statistics_path)
    else:
        print 'No chromosomes to normalise!'
