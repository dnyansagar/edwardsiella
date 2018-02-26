#!/usr/local/python
import os
import re
import sys
import subprocess
import itertools
import argparse
import collections
from Bio import SeqIO
from Bio import pairwise2
from pprint import pprint

def lastz_align(s,q):
	with open('SUBJECT','w') as sf:
		with open('QUERY','w') as qf:
			SeqIO.write(s,sf,'fasta')
			SeqIO.write(q,qf,'fasta')
			sf.close()
			qf.close()
	lastz = subprocess.Popen(['lastz', '--nogapped','--format=general',
			'--nochain','--seed=match10','SUBJECT','QUERY'],
			stdout=subprocess.PIPE)
	score = 0
	for line in lastz.stdout:
		if line.startswith('#'):continue
		elif not re.split('\s+',line):continue
		else:score = int( re.split('\s+',line)[0])
	lastz.stdout.close()
	return score

class PairedSequence(object):
	def __init__(self,seqrec,species,name):
		self.seqrec = seqrec
		self.species = species
		self.pairs = []
	def add_pair(self,paired_sequence):
		self.pairs.append(paired_sequence)
	def __len__(self):
		return len(self.seqrec)
	@property
	def name(self):
		return self.seqrec.name

class SequenceSets(object):
	def __init__(self,seq_dicts):
		self.seq_dicts = seq_dicts
	def find_seq(self,seqid):
		dicts = [sd for sd in seq_dicts if seqid in sd]
		if len(dicts) > 1:
			raise ValueError, 'duplicate id ' + seqid
		if len(dicts) < 1:
			raise ValueError, 'id not found ' + seqid
		return dicts[0][seqid]
		
class Group(object):
	def __init__(self,sequences):
		self.sequences = sequences
		self.length = len(list(sorted(sequences,key=len))[0])
	def __len__(self):
		return self.length

class GroupBuilder(object):
	seq_sets = None
	num_species = 1
	def __init__(self, seq_ids):
		# self.seqs :: Map Species -> List Seq
		self.seqs = collections.defaultdict(list)
		self.representative_seqs = None
		for seq_id in seq_ids:
			seq = seq_sets.find_seq(seq_id)
			self.seqs[seq.species].append(seq)
		self.represented_species = self.seqs.keys()
		# self.seqs is species_name -> PairedSequence objects
		# check that all species are represented
		if len(self.represented_species) < self.num_species:
			return
		# find longest sequence for each species
		map(lambda s: self._build_group_for_species(s), 
			self.represented_species)
		# is self.seqs is list of SeqRecord objects of all species?
		# Where to get pairs? Via PairedSequence or pairs file
	def _build_group_for_species(self,species):
			#look for the longest sequence of species in self.seqs
		longest_seq = list(sorted(self.seqs[species],key=len,reverse=True))[0]
		# look at its pairs, make sure all species are represented
		if len(set(s.species for s in longest_seq.pairs)) < \
			self.num_species - 1:
			return
		# build candidate sequences
		non_paralogs = [s for s in longest_seq.pairs if s.species != species]
		pair_dict = {}
		candidates = None
		for s in non_paralogs:
			if s.species not in pair_dict or \
					len(pair_dict[s.species]) < len(s):
				pair_dict[s.species] = s
		candidates = pair_dict.values()
		# check if the sequence length is sufficent in the pairs found
		seq_covered = lambda s,other: len(other) >= len(s) * .9
		if any(not seq_covered(longest_seq,s) for s in candidates):
			return

		if len(candidates) != self.num_species - 1:
			print len(candidates),self.num_species - 1

		# populate representative_seqs (a list of paired sequences):
		candidate_group = Group(candidates + [longest_seq])
		if self.representative_seqs is not None and \
				len(candidate_group) <= len(self.representative_seqs):
			return

		# run lastz for all in other species
		for sr1,sr2 in itertools.combinations(candidate_group.sequences,2):
			score = lastz_align(sr1.seqrec,sr2.seqrec)
			if score == 0:
				return
		self.representative_seqs = candidate_group

	def get_group(self):
		return self.representative_seqs

''' Usage: script groups pairs cds_seq ...'''

parser = argparse.ArgumentParser(description='Select representative sequences for each species from groupfile obtained via OrthoMCL run')
parser.add_argument('-g','--group', help='group.txt file obtained via OrthoMCL', required=True)
parser.add_argument('-p','--pairs', help='ortholog.txt file in pairs folder of orthomcl run', required=True)
parser.add_argument('-c','--cds', help='path of files with cds sequences', required=True)
args = parser.parse_args()

groups = open(args.group,'r').readlines()
pairs_file = args.pairs
p = args.cds
#len(os.listdir(p))
#GroupBuilder.num_species = 4
GroupBuilder.num_species = int(len(os.listdir(p)))
#groups	= open('/export/home/dnyansagar/Tools/ORTHO/Edwardsiidae/groups_new.txt','r').readlines()
#pairs_file	= '/export/home/dnyansagar/Tools/ORTHO/Edwardsiidae/orthomcl/pairs/orthologs.txt'
#p = '/home/dnyansagar/Tools/ORTHO/Edwardsiidae/files/cds/'

def make_species_dict(species_name):
	return {sr.name: PairedSequence(sr,species_name,sr.name)
			for sr in SeqIO.parse(os.path.join(p,'{}.fasta'.format(species_name)),
								  'fasta')}

par_seqs = make_species_dict('Par') 
eli_seqs = make_species_dict('Eli') 
eca_seqs = make_species_dict('Eca') 
nve_seqs = make_species_dict('Nve') 

seq_dicts = [eca_seqs,eli_seqs,nve_seqs,par_seqs]

# create the sequence sets - this object indexes all sequences 
# and groups them by species (the "sets")
seq_sets = SequenceSets([eca_seqs,eli_seqs,nve_seqs,par_seqs])

# using the pairs file, look for sequences in the seqeunce sets
# and add the information as to which sequences pair off to each other
for line in open(pairs_file):
	id1,id2 = line.split()[:2]
	seq1 = seq_sets.find_seq(id1)
	seq2 = seq_sets.find_seq(id2)
	seq1.add_pair(seq2)
	seq2.add_pair(seq1)

for line in groups:
	spl = re.split('\s+',line.strip())
	gb = GroupBuilder(spl[1:])
	group = gb.get_group()
	if group is not None:
		print '\t'.join(s.name for s in group.sequences)
		
