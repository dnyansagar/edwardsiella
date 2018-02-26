#!/usr/local/python
import os,re,sys,subprocess
from Bio import SeqIO
from Bio import pairwise2
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
match = open('/export/home/dnyansagar/Tools/ORTHO/Edwardsiidae/pair_list.txt','r').readlines()
#match = open( '/export/home/dnyansagar/Tools/ORTHO/Edwardsiidae/Filtered_groups_15Apr_reverse.txt','r').readlines()
p = '/home/dnyansagar/Tools/ORTHO/Edwardsiidae/files/cds/'

def getfile(path,string):
        for fname in os.listdir(path):
                if fname.endswith('.fasta'):
                        if string in open(path+fname).read():
                                return fname
                        else:
                                continue
def getseq(description,filename):
        fasta=SeqIO.parse(p+filename,"fasta")
        for record in fasta:
                if description in record.description:
                        return record
                else:
                        continue

def lastz(file1,file2):
        lastz = subprocess.Popen(['lastz', '--format=general', '--chain','--noytrim', 
								file1,file2],stdout=subprocess.PIPE)
        for line in lastz.stdout:
                if line.startswith('#'):continue
                elif not re.split('\s+',line):continue
                else:
                        spl2 = re.split('\s+',line)
                        ide = spl2[12].replace('%','')
                        return float(ide)

p1 = list();p2 = list();p3 = list();p4 = list();p5 = list();p6 = list()
for line in match:
        spl = re.split('\s+',line)
	ecaid = spl[1].split('|')[0].replace("'",""); 
	eliid = spl[2].replace("'","");
	nveid = spl[3].replace("'","");
	parid = spl[4].split('|')[0].replace("'","")
	fl1 = getfile(p,ecaid);fl2 = getfile(p,eliid);fl3 = getfile(p,nveid);fl4 = getfile(p,parid)
	i1 = getseq(ecaid,fl1);j1 = getseq(eliid,fl2);k1 = getseq(nveid,fl3);l1 = getseq(parid,fl4); 
        SeqIO.write(i1,'ec_fasta',"fasta");SeqIO.write(j1,'el_fasta',"fasta");SeqIO.write(k1,'nv_fasta',"fasta");SeqIO.write(l1,'pa_fasta',"fasta"); 
	
	lastz1 = lastz('ec_fasta','el_fasta')
	if lastz1 is not None:
		print lastz1
		p1.append(lastz1)
	
	lastz2 = lastz('ec_fasta','pa_fasta')
	if lastz2 is not None:
		print lastz2
		p2.append(lastz2)
	
	lastz3 = lastz('pa_fasta','nv_fasta')
	if lastz3 is not None:
		print lastz3
		p3.append(lastz3)
	
	lastz4 = lastz('nv_fasta','el_fasta')
	if lastz4 is not None:
		print lastz4
		p4.append(lastz4)
		
	lastz5 = lastz('pa_fasta','el_fasta')
	if lastz5 is not None:
		print lastz5
		p5.append(lastz5)
	
	lastz6 = lastz('nv_fasta','ec_fasta')
	if lastz6 is not None:
		print lastz6
		p6.append(lastz6)
	
print len(p1),float(sum(p1)/len(p1))
print len(p2),float(sum(p2)/len(p2))
print len(p3),float(sum(p3)/len(p3))
print len(p4),float(sum(p4)/len(p4))
print len(p5),float(sum(p5)/len(p5))
print len(p6),float(sum(p6)/len(p6))
