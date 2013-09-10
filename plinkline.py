#!/usr/bin/python
#THE PIPELINE! - things to note! Everywhere this assumes you have 39 chromosomes. Replacing all the range(1,40) and range(39) lines with range(1,N+1) or range(N) for N chromosomes should work, but you may need to create some new directories to house them.
# Welcome to the program! This file basically provides the framework for the pipeline, calling functions from other files (which should be in the same directory!).


# Importing some basics

import os
from multiprocessing import Pool
import glob
import sys
import argparse
from time import localtime, strftime
print 'The time is', strftime("%a, %d %b %Y %X +0000", localtime())

import locations                # locations.py contains the SYSTEM SPECIFIC directories for the various inputs and outputs - you will need to edit this if you run this pipeline on a different computer, sorry!


# Defining file locations - as above, we're just fetching these from locations.py.

ped_folder_path = locations.ped_folder_path                     # location of ped files
master_ped_name = locations.master_ped_name                     # prefix of .chrN.ped files containing the data - for example, for my analyses it was 'phase2+other+pfizer31-36.cf3.1'
keeplists_folder_path = locations.keeplists_folder_path         # this folder will contain files whose names are input to the program, containing individual IDs (found in ped file)
ihsmap_folder_path = locations.ihsmap_folder_path               # this folder contains the ihsmaps
ihshap_folder_path = locations.ihshap_folder_path               # folder containing ihshaps (output of translation algo)
statistics_folder_path = locations.statistics_folder_path       # folder containing output of statistical tests
regions_folder_path = locations.regions_folder_path             # folder containing regions of interest - this is where the final product will end up!

ihs_path = locations.ihs_path                                   # path of ihs script (these are the executables from the Pritchard software)
xpehh_path = locations.xpehh_path                               # path of xpehh script

# Some more constants

snps_per_bin = '25'                                             # this is how many SNPs to put in a window for the purpose of region-finding... this is a string because it gets passed to Rpy2 later

# Parse command line options

parser = argparse.ArgumentParser(description='Welcome to the limit!')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--ihs', nargs=1, help='run integrated haplotype score (iHS) test on a single population', metavar='filename')
group.add_argument('--xpehh', nargs=2, help='run cross-population extended haplotype homozygosity (XP-EHH) test between two populations',metavar=('filename1','filename2'))
group.add_argument('--extract_list', nargs=1, help='batch extract a list of populations', metavar='list_filename')

parser.add_argument('--extract', help='perform extraction component (ped -> breed_ped)', action='store_true')
parser.add_argument('--translate', help='perform translation component (breed_ped -> ihshap)',action='store_true')
parser.add_argument('--run_test', help='perform test componoent (ihshap+ihsmap -> test_stats',action='store_true')
parser.add_argument('--normalise', help='perform normalisation component (test_stats -> test_stats_normed)',action='store_true')
parser.add_argument('--analyse', help='perform analysis component (test_stats_normed -> regions)',action='store_true')
parser.add_argument('--chro','--chr', help='specify to a single chromosome. note! some options will always require all chromosomes (normalisation, analysis), so this option does not affect them')
parser.add_argument('--cores', help='how many cores to use? (default is 4)', choices=['4','8'], default = 4)

args = parser.parse_args()

# Translate into program variables

options = [key for key in vars(args) if type(vars(args)[key])== type([]) or vars(args)[key] not in {None, False}]         # this is a bit of a mess, but it seems to work - note! seems not to work in some versions of python, or something... this can just be commented out if it seems to be causing problems!

extract_flag = args.extract
translate_flag = args.translate
run_test_flag = args.run_test
normalise_flag = args.normalise
analyse_flag = args.analyse

# Specify a certain chromosome - note! As above, some options always require all chromosomes (normalisation, analysis), so this option won't affect them.
# Here is a more elegant way to generate these lists:
# chr_lists = [ [chr_list[x] for x in range(y,39,cores)] for y in range(0,cores)]
# And then chr_list1 = chr_lists[0], and so on and so forth
# The seemingly-random ordering of the chromosomes in the lists here (and the reason I didn't use the more elegant method above) is that they're ordered by size, so that no process should finish significantly faster than any of the rest. -> Mild optimisation.

chr_opt = args.chro
if chr_opt == None:
    cores = int(args.cores)
    if cores == 8:
#        print 'Using 8 cores!'
        chr_list1 = [36,29,23,11,12]
        chr_list2 = [32,34,14,9]
        chr_list3 = [37,30,18,15,2]
        chr_list4 = [38,31,21,22,4]
        chr_list5 = [33,19,20,17,5]
        chr_list6 = [28,39,10,6,3]
        chr_list7 = [26,27,25,13,7]
        chr_list8 = [35,24,16,8,1]
    else:
#        print 'Using 4 cores!'
        chr_list1 = [36,33,29,19,23,20,11,17,12,5]
        chr_list2 = [38,32,28,34,39,14,10,9,6,3]
        chr_list3 = [37,26,30,27,18,25,15,13,2,7]
        chr_list4 = [35,31,24,21,16,22,8,4,1]
        chr_list5 = []
        chr_list6 = []
        chr_list7 = []
        chr_list8 = []
else:
    chrom = int(args.chro)
    chr_list1 = [chrom]
    chr_list2 = []
    chr_list3 = []
    chr_list4 = []
    chr_list5 = []
    chr_list6 = []
    chr_list7 = []
    chr_list8 = []

# Which test?

if args.xpehh!=None:
    test_type = 'xpehh'

    pop_1_name = args.xpehh[0]
    pop_2_name = args.xpehh[1]
    pop_names = [pop_1_name,pop_2_name]

    test_data_name = pop_1_name+'-'+pop_2_name

elif args.ihs!=None:
    test_type = 'ihs'

    pop_name = args.ihs[0]
    pop_names = [pop_name]

    test_data_name = pop_name

else:
    extract_flag = True
    pop_names =[]
    extraction_list = args.extract_list[0]
    for line in open(keeplists_folder_path+extraction_list,'r'):
        pop_names.append(line.strip())

# Parallelise!

p = Pool(cores)                                     # Cores is either 4 or 8. Default is 4! To make this all magically work on more cores, the chr_lists need to be changed. See above (where I defined the chr_lists for the first time) for a way to do this.

chr_lists_all = [chr_list1,chr_list2,chr_list3,chr_list4,chr_list5,chr_list6,chr_list7,chr_list8]
chr_lists = chr_lists_all[:cores]                   # As above, if you've picked 4 cores, the last 4 chr_lists are just empty, so we cut them out.

# Importing functions from elsewhere! - if these files are not in the same folder as this script, trouble is ahead.

import extract
import translate
import run_test
import norm
import analyse

# The main routine

optional_flags = {extract_flag,translate_flag,run_test_flag,normalise_flag,analyse_flag}

if True in optional_flags:
    print 'Running with options:', ', '.join(options)+' '+str(cores)
else:                                           # This is how no-flags implies all-flags.
    extract_flag = True
    translate_flag = True
    run_test_flag = True
    normalise_flag = True
    analyse_flag = True
    print 'Running full analysis with test:', test_type

# ---- create per-population ped files ---- #
# This just extracts lines from a master ped file which have IDs in the supplied lists, and puts them into their own ped file. Nothing too fancy.
if extract_flag ==True:
    extract.run_extraction(pop_names, keeplists_folder_path, chr_lists, ped_folder_path, master_ped_name, p,cores)

# ---- convert ped file into ihshap file ---- #
# See the Pritchard scripts documentation for details about the format of ihshap files. Essentially just the genotype information from ped files, but base pairs recoded to 0 or 1 depending on whether or not they're ancestral. In addition, we cut out any SNPs for which our reference (culpeo fox) was not homozygous.
if translate_flag == True:
    translate.run_translation(chr_lists,pop_names,ped_folder_path,ihshap_folder_path,ihsmap_folder_path,p,cores)

# ---- calculate the statistic! ---- #
# Running the Pritchard scripts.
if run_test_flag == True:
    run_test.run_the_test(pop_names, chr_lists, test_type, statistics_folder_path, ihshap_folder_path, ihsmap_folder_path, ihs_path, xpehh_path, p,cores)

# ---- normalise the test statistic ---- #
# The output of the Pritchard scripts are non-normalised, so here we turn them to standard normals. iHS data gets normalised in 5% frequency bins, both are genome-wide.
if normalise_flag == True:
    norm.run_normalisation(test_data_name,test_type,statistics_folder_path)

# ---- using the test statistic, find putative regions of selection ---- #
# We use sliding windows across each chromosome. Each window is given a score, and we take the top 1% of the distribution of this score (across all chromosomes) to decide on candidate windows.
#   The output is a file whose columns are:
# chr   region_start    region_end  region_statistic    top_SNP location_of_top_SNP test_statistic_for_top_SNP
if analyse_flag == True:
    analyse.find_regions(test_data_name, test_type, snps_per_bin, statistics_folder_path, regions_folder_path)

# ---- create a logfile ---- #

