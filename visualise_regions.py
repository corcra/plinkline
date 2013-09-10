#!/usr/bin/python
# This script generates an r script which can be run to look at a specific chromosome of a dataset, highlighting regions which have been identified as possible targets of selection.
# Therefore it requires: normed test statistic data, and a .regions file.
# It sticks the R script in test_statistics/data-type (eg xpehh or ihs), from which you can run it in R with the command source('look_chr.r').

import sys
import locations

#argv[1] dataset
#argv[2] data type
#argv[3] chr

if len(sys.argv)<3:
    sys.exit('dataset, datatype, chr')

job_name = sys.argv[1]
data_type = sys.argv[2]

if data_type == 'ihs':
    col = '11'
elif data_type == 'xpehh':
    col = '6'
else:
    sys.exit('Datatype (second argument) must be xpehh or ihs')

chro = int(sys.argv[3])

#acquiring the regions
regions_file = locations.regions_folder_path+job_name+'.regions'
regions = []
for line in open(regions_file,'r'):
    chrom = int(line.split()[0])
    if chrom == chro:
        lower = line.split()[1]
        upper = line.split()[2]
        frac = line.split()[3]
        regions.append((lower,upper,frac))


#setting up files
basic_path = locations.statistics_folder_path+data_type+'/'
filename = basic_path+'chr'+str(chro)+'/'+job_name+'.chr'+str(chro)+'.norm'
outfile_path = locations.statistics_folder_path+data_type+'/look_chr.r'
outfile = open(outfile_path,'w')


outfile.write('data<-read.table(\''+filename+'\',as.is=TRUE)\n')
outfile.write('data<-data[is.finite(data[,'+col+']),]\n')

outfile.write('library(plotrix)\nlibrary(gplots)\n')
outfile.write('par(family="Avenir",bty="l")\n')
outfile.write('x<-seq(0,max(data[,2]),max(data[,2])/1000)\nx<-x[1:1000]\ny<-rep(c(-2.5,2.5),500)\n')
outfile.write('ymax<-max(abs(data[,'+col+']),na.rm=TRUE)\n')

if data_type == 'ihs':
    outfile.write('plot(x/(1e07),y,type="p",pch=".",col="darkorange1",ylim=c(0,1.03*ymax),main=paste("'+job_name+': chr '+str(chro)+'"),xlab="Physical Distance (e+07)",ylab="|iHS|")\n')
    absornot = 'abs'

if data_type == 'xpehh':
    outfile.write('plot(x/(1e07),y,type="p",pch=".",col="darkorange1",ylim=c(1.03*min(data[,'+col+'],na.rm=TRUE),1.03*max(data[,'+col+'],na.rm=TRUE)),main=paste("'+job_name+': chr '+str(chro)+'"),xlab="Physical Distance (e+07)",ylab="XP-EHH")\n')
    absornot = ''

start = 1
for region in regions:
    lower = region[0]
    upper = region[1]

    #first the non-region up to this region
    end = lower
    outfile.write('region<-data[(data[,2]>='+str(start)+')&(data[,2]<'+str(end)+'),]\n')
    outfile.write('points(region[,2]/(1e07),'+absornot+'(region[,'+col+']),col="black",type="h",lw=0.8)\n')

    #next the actual region
    start = lower
    end = upper
    outfile.write('region<-data[(data[,2]>='+str(start)+')&(data[,2]<'+str(end)+'),]\n')
    outfile.write('points(region[,2]/(1e07),'+absornot+'(region[,'+col+']),col="darkorange1",type="h",lw=0.8)\n')

    #now prepare it for the next region
    start = str(int(end)+1)

outfile.write('region<-data[(data[,2]>='+str(start)+'),]\n')
outfile.write('points(region[,2]/(1e07),'+absornot+'(region[,'+col+']),col="black",type="h",lw=0.8)\n')

#outfile.write('legend(0,1.05*ymax,"Green : top 1% of windows wrt fraction of SNPS with |'+data_type+'|>2.5",bg="mediumpurple1",cex=0.75)\n')
