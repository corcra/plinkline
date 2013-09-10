# This is all about finding regions after the iHS/XP-EHH statistics have been calculated.
# Here is how the region-finding algorithm works:
#   Split the chromsome into 25-SNP bins. 
#   Do sliding windows around each of these (so the first bin is SNPs 1-25, and the second window is SNPs 2-26, etc. This system is why there is so much machinery for dealing with edge-cases at the end of the chromosome.)
#   For each WINDOW a score is assigned. This is either 'proportion of extreme SNPs', where extreme means a test statistic >2.5 (abs value), or 'average SNP value'. At the time of writing, the latter is in force.
#   For each BIN, the window with the greatest score is taken. In addition, the most extreme SNP in that region is also recorded.
#   We then have a sample of maximal-regions, one for each bin. 
#   This is repeated for all chromosomes.
#   An empirical (genome_wide) distribution is created over these, and the top 1% taken as putative regions of selection. Thus the algorthim will always output at least (num_bins)*(num_chr)/100 regions, even if some of these are not very interesting. (eg if the greatest SNP score across the whole genome (data-permitting) were 2, this is not very 'unusual' for a normally-distributed value. However, in the absence of some method to empirically produce a distribution function for these scores, I would be hesitant to try assign any p-values.)

# NOTE!! This script uses rpy2, which is an R-interface for python. rpy2 requires R version 2.16 or newer to be installed. The cac computers have R 2.10 installed. This is a problem. I have tried rectifying this by locally installing R 3.0.1 and reinstalling rpy2 (found in my home directory) while directing it to the local R installation, but this is *inexplicably not working*. Short of rewriting this entire script in numpy or straight-up R, I am out of ideas.

import sys
from subprocess import call
import rpy2.robjects as robj
import os
import numpy as np

def overlap(nuple_1,nuple_2):
    lower_1 = int(nuple_1[0])
    lower_2 = int(nuple_2[0])

    upper_1 = int(nuple_1[1])
    upper_2 = int(nuple_2[1])

    if lower_1<=lower_2 and lower_2<=upper_1:
        return True
    elif lower_2<=lower_1 and lower_1<=upper_2:
        return True
    else:
        return False

def get_region_info(data,index,snps_per_bin,chro,i):
    bin_data = data[(index+i):(index+i+int(snps_per_bin))]               # this is ever so slightly different to the R way
    lower = bin_data[0][1]
    upper = bin_data[int(snps_per_bin)-1][1]
    snp_scores = [snp_data[2] for snp_data in bin_data]
#    snp_scores = [bin_data[n][2] for n in range(len(bin_data))]
#   extreme_snps = [abs(snp_score)>=2.5 for snp_score in snp_scores]   # the 'def' method
#   prop = np.mean(extreme_snps)
    prop = np.mean(map(abs,snp_scores))                                 # the 'abs' method
    max_snp = sorted(bin_data, key=lambda snp_data: abs(snp_data[2]),reverse=True)[0]
    max_snp_id = max_snp[0]
    max_snp_loc = max_snp[1]
    max_snp_stat = max_snp[2]
    region_info = [int(chro), lower, upper, prop, max_snp_id, max_snp_loc, max_snp_stat]
    return region_info

def numpy_part(test_data_name, test_type, col, snps_per_bin, statistics_folder_path):       #it really doesn't use numpy, so this naming is deceptive

    regions = []
    for chro in range(1,40):
        chr_data = []
        test_stat_file_path = statistics_folder_path+test_type+'/chr'+str(chro)+'/'+test_data_name+'.chr'+str(chro)+'.norm'
        for line in open(test_stat_file_path,'r'):
            snpid = line.split()[0]
            loc = int(line.split()[1])
            test_statistic = line.split()[int(col)-1]
            if not test_statistic in {'inf', 'nan'}:
                snp_data = (snpid, loc, float(test_statistic))
                chr_data.append(snp_data)
        print 'Chr'+str(chro)+' read in.'
        
        length = len(chr_data)
        if length<=int(snps_per_bin):
            print 'Not enough data on chromosome '+str(chro)+' ('+str(length)+') - skipping!'
            continue

        num_bins = int(np.floor(length/float(snps_per_bin)))

        print 'Analysing chr'+str(chro)+'...'
        index = 0
        for j in range(num_bins-1):
            index = j*int(snps_per_bin)
            local_region = []
            for i in range(int(snps_per_bin)):
                region_info = get_region_info(chr_data,index,snps_per_bin,chro,i)
                local_region.append(region_info)
            max_local_region = sorted(local_region, key = lambda region: region[3], reverse=True)[0]
            regions.append(max_local_region)
            
        print 'Finishing off chr'+str(chro)
        local_region = []
        final_index = (num_bins-1)*int(snps_per_bin)
        leftover = length - final_index
        maxrange = leftover - int(snps_per_bin)
        for i in range(maxrange+1):
            region_info = get_region_info(chr_data,final_index,snps_per_bin,chro,i)
            local_region.append(region_info)
        max_local_region = sorted(local_region, key = lambda region: region[3])[0]
        regions.append(max_local_region)

    means = [region[3] for region in regions]
    num_regions = len(regions)
    top_1_percent = int(np.floor(num_regions/100.0))
    top_regions = sorted(regions, key = lambda region: region[3],reverse=True)[0:top_1_percent]
    cutoff = top_regions[-1][3]
    print 'Cutoff: ', cutoff

    print 'Giving it back...'                         # I could do this with a list comprehension, but i want the nice names

    chro_l = []
    lower_l = []
    upper_l = []
    prop_l = []
    max_stat_id_l = []
    max_stat_loc_l = []
    max_stat_stat_l = []

    for region in top_regions:
        chro_l.append(region[0])
        lower_l.append(region[1])
        upper_l.append(region[2])
        prop_l.append(region[3])
        max_stat_id_l.append(region[4])
        max_stat_loc_l.append(region[5])
        max_stat_stat_l.append(region[6])

    lists = [chro_l, lower_l, upper_l, prop_l, max_stat_id_l, max_stat_loc_l, max_stat_stat_l]
    return lists


def r_part(test_data_name, test_type, col, snps_per_bin, statistics_folder_path):           # since i wrote the numpy part, this section is semi-defunct, but it may be useful to compare, in case of slight differences (which exist)

	#this is all quite horrible and inefficient
    robj.r('regions<-cbind(NA,NA,NA,NA,NA,NA,NA)')
    for chro in range(1,40):
        print 'Reading chr'+str(chro)+' into R(py2).'
        test_stat_file_path = statistics_folder_path+test_type+'/chr'+str(chro)+'/'+test_data_name+'.chr'+str(chro)+'.norm'

        dat= robj.r('data<-read.table("'+test_stat_file_path+'",as.is=TRUE)')               # note, this does not work well with empty files, but hopefully it won't come to that
        robj.r('data<-data[is.finite(data[,'+col+']),]')

        ell = robj.r('l<-length(data[,2])')
        if ell[0]<=int(snps_per_bin):
            print 'Not enough data on chromosome '+str(chro)+' ('+str(ell[0])+') - skipping!'
            continue

        robj.r('chro<-'+str(chro))
        num = robj.r('num_bins<-floor(l/('+snps_per_bin+'))')
        num_bins = int(num[0])
        print 'Analysing chr'+str(chro)+'...'
        robj.r('index<-0')
        for j in range(num_bins-1):
            robj.r('index<-'+str(j)+'*'+snps_per_bin)
            robj.r('local_region<-cbind(NA,NA,NA,NA,NA,NA,NA)')
            for i in range(1,int(snps_per_bin)+1):
                robj.r('bin<-data[(index+'+str(i)+'):(index+'+str(i)+'+'+snps_per_bin+'),]')
                robj.r('lower<-bin[1,2]')
                robj.r('upper<-bin['+snps_per_bin+',2]')

# ---- here we choose the method of scoring windows ----- #
#                robj.r('prop<-mean(abs(bin[,'+col+'])>2.5)') # 'default'
                robj.r('prop<-mean(abs(bin[,'+col+']),na.rm=TRUE)')         # abs'


            #SNP with max test statistic
                robj.r('max_stat<-bin[sort.list(abs(bin[,'+col+']),decreasing=TRUE),][1,]')
            #this is clearly a mess
                robj.r('max_stat_id<-max_stat[1,1]')
                robj.r('max_stat_loc<-max_stat[1,2]')
                robj.r('max_stat_stat<-max_stat[1,'+col+']')
                robj.r('region_info<-cbind(chro,lower,upper,prop,max_stat_id,max_stat_loc,max_stat_stat)')
                robj.r('local_region<-rbind(local_region,region_info)')
        #choosing max region by prop score
            robj.r('local_region<-local_region[2:length(local_region[,1]),]')
            robj.r('max_local_region<-local_region[sort.list(local_region[,4],decreasing=TRUE),][1,]')
            robj.r('regions<-rbind(regions,max_local_region)')

        print 'Finishing off chr'+str(chro)
        robj.r('local_region<-cbind(NA,NA,NA,NA,NA,NA,NA)')
        robj.r('final_index<-(num_bins-1)*'+snps_per_bin)
        leftover = robj.r('leftover<-(l-final_index)')
        mrange = robj.r('leftover-'+snps_per_bin)
        maxrange = int(mrange[0])
        for i in range(1,maxrange+2):
            robj.r('bin<-data[(final_index+'+str(i)+'):(final_index+'+str(i)+'+'+snps_per_bin+'),]')
            robj.r('lower<-bin[1,2]')
            robj.r('upper<-bin['+snps_per_bin+',2]')
#            robj.r('prop<-mean(abs(bin[,'+col+'])>2.5)')        # 'default'
            robj.r('prop<-mean(abs(bin[,'+col+']),na.rm=TRUE)') # 'abs'
            robj.r('max_stat<-bin[sort.list(abs(bin[,'+col+']),decreasing=TRUE),][1,]')
        #this is clearly a mess
            robj.r('max_stat_id<-max_stat[1,1]')
            robj.r('max_stat_loc<-max_stat[1,2]')
            robj.r('max_stat_stat<-max_stat[1,'+col+']')
            robj.r('region_info<-cbind(chro,lower,upper,prop,max_stat_id,max_stat_loc,max_stat_stat)')
            robj.r('local_region<-rbind(local_region,region_info)')
        # the array was initialised with a row of NA so we're getting rid of those now. If the array only has two rows, we take the second as the max region (it is the only region).
        robj.r('if (length(local_region[,1])>2){local_region<-local_region[2:length(local_region[,1]),]\nmax_local_region<-local_region[order(local_region[,4],decreasing=TRUE),][1,]} else{max_local_region<-local_region[2,]}')
        robj.r('regions<-rbind(regions,max_local_region)')
    #cut off the initial list of NAs (both local_region and regions had initial list of NAs... there is probably a way to initialise a dataframe that doesn't require this, but I don't know it... yet
        robj.r('regions<-regions[2:length(regions[,1]),]')

    robj.r('means<-regions[,4]')
    robj.r('CDF<-ecdf(means)')
    robj.r('cutoff<-quantile(CDF,0.99)')
    cutoff = robj.r('cutoff')
    print 'Cutoff: ',cutoff[0]

    robj.r('good_regions<-regions[regions[,4]>=cutoff,]')

    print 'Giving it back to python...'                         # I could do this with a list comprehension, but i want the nice names

    chro_l = robj.r('good_regions[,1]')
    lower_l = robj.r('good_regions[,2]')
    upper_l = robj.r('good_regions[,3]')
    prop_l = robj.r('good_regions[,4]')
    max_stat_id_l = robj.r('good_regions[,5]')
    max_stat_loc_l = robj.r('good_regions[,6]')
    max_stat_stat_l = robj.r('good_regions[,7]')

    lists = [chro_l, lower_l, upper_l, prop_l, max_stat_id_l, max_stat_loc_l, max_stat_stat_l]
    return lists

def py_part(lists):
    print 'Thanks R, python is go!'

    chro_l = lists[0]
    lower_l = lists[1]
    upper_l = lists[2]
    prop_l = lists[3]
    max_stat_id_l = lists[4]
    max_stat_loc_l = lists[5]
    max_stat_stat_l = lists[6]

    master_list = []
    for i in range(1,40):
        master_list.append([])

    length = len(chro_l)
    #populate list (separating chromosomes)
    for i in range(length):
        if str(chro_l[i]) == '0' or str(chro_l[i]) == 'NA':
            print 'Problematic chromosome number detected (chr'+str(chro_l[i])+') - skipping!', i
            continue
        chro = int(chro_l[i])
        lower = int(lower_l[i])
        upper = int(upper_l[i])
        prop = float(prop_l[i])
        max_stat_id = str(max_stat_id_l[i])
        max_stat_loc = int(max_stat_loc_l[i])
        max_stat_stat = float(max_stat_stat_l[i])
        master_list[chro-1].append((lower,upper,prop,max_stat_id,max_stat_loc,max_stat_stat))

    print 'Master list created.'
    regions_to_keep = dict()
    print 'Analysing...'
    for i in range(1,40):
        regions = set()
        for nuple in master_list[i-1]: #the ith chromosome is the ith-1 element in master_list, remember!
            rem = set()
            add = 0
            unique = 1
            this_stat = nuple[2]
            for region in regions:
                if overlap(nuple,region):
                    unique = 0
                    region_stat = region[2]
                    if this_stat > region_stat:
                        add = 1
                        rem.add(region)
            if unique == 1 or add == 1:
                regions.add(nuple)
            regions_temp = regions.difference(rem)
            regions = regions_temp
        if len(regions)!=0:
            regions_to_keep[i] = sorted(regions) #if the lower region is the first element of the tuple, this should sort our regions nicely

    return regions_to_keep

def save_regions(regions_folder_path, test_data_name, regions_to_keep):
    out_file_path = regions_folder_path+test_data_name+'.regions'
    out_file = open(out_file_path,'w')
    print 'Recording regions of putative selection to '+out_file_path
    for chro in sorted(regions_to_keep):
        for region in regions_to_keep[chro]:
            lower = str(region[0])
            upper = str(region[1])
            prop = str(region[2])
            max_stat_id = str(region[3])
            max_stat_loc = str(region[4])
            max_stat_stat = str(region[5])
            out_file.write(str(chro)+' '+lower+' '+upper+' '+prop+' '+max_stat_id+' '+max_stat_loc+' '+max_stat_stat+'\n')

def find_regions(test_data_name, test_type, snps_per_bin, statistics_folder_path, regions_folder_path):             # this is not very modularised, clearly
    if test_type == 'ihs':
        col='11'
    elif test_type == 'xpehh':
        col='6'

    if not os.path.exists(statistics_folder_path+test_type+'/chr1/'+test_data_name+'.chr1.norm'):
    	print 'Please normalise test statistic first! - cannot look for regions!'
    else:
        print 'Looking for regions of selection in '+test_data_name
#    	lists = r_part(test_data_name, test_type, col, snps_per_bin, statistics_folder_path)                # using R
    	lists = numpy_part(test_data_name, test_type, col, snps_per_bin, statistics_folder_path)            # using numpy (mostly python)
    	regions_to_keep = py_part(lists)
    	save_regions(regions_folder_path, test_data_name, regions_to_keep)
