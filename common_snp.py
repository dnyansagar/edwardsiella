#!/usr/local env python
import re,sys,os
import argparse
parser = argparse.ArgumentParser(description='Compare SNP locations in varscan filtered output')
parser.add_argument('-ort','--ortholog', help='Ortholog file from OrthoMCl', required=True)
parser.add_argument('-eca','--eca_snp', help='snp file for E. carnea', required=True)
parser.add_argument('-eli','--eli_snp', help='snp file for E. lineata', required=True)
parser.add_argument('-par','--par_snp', help='snp file for E. parasite', required=True)
args = parser.parse_args()

""" Files Used
Orthologs ='/export/home/dnyansagar/Tools/ORTHO/Edwardsiidae/Filtered_groups_15Apr.txt'
vcf_eca = '/export/home/dnyansagar/1/SNP_LOCUS/varscan/Eca_0.01.filtered.varscan'
vcf_eli = '/export/home/dnyansagar/1/SNP_LOCUS/varscan/Eli_0.01.filtered.varscan'
vcf_par = '/export/home/dnyansagar/1/SNP_LOCUS/varscan/Par_0.01.filtered.varscan'
"""
# Get lines from the file containing given text
def getlines(file1,txt):
        lst = list()
        for line in open(file1,'r'):
                if txt in line:
                        lst.append(line)
        return lst #''.join(lst)

eca_eli=list();	eca_par=list();	eli_par=list();	eca_eli_par=list()

def compareSNP(ecline, elline, paline):
	ecsnp = re.split('\s+',ecline);	elsnp = re.split('\s+',elline);	pasnp = re.split('\s+',paline)
	if  ecsnp[1] == elsnp[1]:
		mut = ecsnp[0]+'-'+ecsnp[1]+'_'+elsnp[0]+'-'+elsnp[1]
		if mut not in eca_eli: eca_eli.append(mut)
	if ecsnp[1] == pasnp[1]:
		mut = ecsnp[0]+'-'+ecsnp[1]+'_'+pasnp[0]+'-'+pasnp[1]
		if mut not in eca_par:eca_par.append(mut)
	if elsnp[1] == pasnp[1]:
		mut  = elsnp[0]+'-'+elsnp[1]+'_'+pasnp[0]+'-'+pasnp[1]
		if mut not in eli_par: eli_par.append(mut)
	if ecsnp[1]==elsnp[1]==pasnp[1]:
		mut = ecsnp[0]+'-'+ecsnp[1]+elsnp[0]+'-'+elsnp[1]+'_'+pasnp[0]+'-'+pasnp[1]
		if mut not in eca_eli_par:eca_eli_par.append(mut)

for x in open(args.ortholog,'r').readlines():
	spl = re.split('\s+',x)
	ec_acc = spl[0].split('|')[1].replace('c','')			# Get E.carnea id and format it to match id in Orthologs
	ec = getlines(args.eca_snp,ec_acc)
	el_acc = spl[1].split('|')[1].replace('lineata_','')	# Get E.lineata id and format it to match id in Orthologs
        el = getlines(args.eli_snp, el_acc)
	pa_acc = spl[3].split('|')[1]							# Get E.parasite id and format it to match id in Orthologs
        pa = getlines(args.par_snp, pa_acc)
	for x, y ,z in [(x,y,z) for x in ec for y in el for z in pa]:
		abcd = compareSNP(x,y,z)
print 'E. carnea - E.lineata:', len(eca_eli),'\n',\
		'E.carnea - E. parasite:',len(eca_par),'\n',\
		'E. lineata - E. parasite',len(eli_par),'\n',\
		'Common to all:', '\t',len(eca_eli_par)
