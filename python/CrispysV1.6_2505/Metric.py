import pickle
import numpy as np
from functools import reduce
import UPGMA
import random
import string
from scipy.spatial.distance import minkowski
import sys
sys.path.append("/bioseq/crista/CRISTA_online/")
from main import check_input_vars, get_features, predict_on_df
import numpy as np

def pos_in_metric_general(t, df, base, cfd_dict = None):
	'''
	:param t: target
	:param df: distance function
	:param base: a list of strings, spreading the space. each string is an sgRNA
	:return: a vector of distances between those strings
	'''
	#return list(map(lambda sg: df(sg,t), base))
	if cfd_dict:
		Vetorize = np.vectorize(lambda sg: df(sg, t, cfd_dict))
	else:
		Vetorize = np.vectorize(lambda sg: df(sg,t))
	return Vetorize(base)

def pos_in_metric_cfd(t, cfd_dict = None):
	'''
	:param t: target
	 implement a version of the cfd score, in which
	:return:
	'''
	if not dicti:
		dicti = pickle.load(open("cfd_dict.p",'rb'))
	#print(dicti)
	Nucs = ['A','C','G', 'U']
	point = [0 for i in range(len(t)*len(Nucs))]
	i=0
	for pos in range(len(t)):
		for Nuc in Nucs:
			key = ('r'+ Nuc +':d'+ t[pos], pos+1)
			if key in dicti:
				point[i] = dicti[('r'+ Nuc +':d'+ t[pos], pos+1)]
			else:
				point[i] = 1
			i += 1
	return point

def pos_in_metric_cfd_np(t, dicti):
	'''
	:param t: target
	 implement a version of the cfd score, in which
	:return:


	there is a bug here - the code and the dictinary dose not fit.
	'''
	if not dicti:
		dicti = pickle.load(open("/groups/itay_mayrose/galhyams/MULTICRISPER/codeV1.3SetCover/cfd_dict.p",'rb'))
		#dicti = pickle.load(open("D:\\Lab\\BackUp\\cfd_dict.p",'rb'))

	Nucs = ['A','C','G', 'U']
	point = np.zeros(len(t)*len(Nucs))
	#point = [0 for i in range(len(t)*len(Nucs))]
	i=0
	for pos in range(len(t)):
		for Nuc in Nucs:
			key = ('r'+ Nuc +':d'+ t[pos], pos+1)
			if key in dicti:
				point[i] = dicti[('r'+ Nuc +':d'+ t[pos], pos+1)]
			else:
				point[i] = 1
			i += 1
	return point
###scoring funct: CFD and CRISTA###
def cfd_funct(sgRNA, target, dicti = None):
	'''my implementation of this function'''
	if not dicti:
	#dicti = pickle.load(open("D:\\Lab\\BackUp\\cfd_dict.p",'rb'))
		dicti = pickle.load(open("/groups/itay_mayrose/galhyams/CrispysV1.6/cfd_dict.p",'rb'))

	#print(('r'+sgRNA[0]+':d'+target[0],1))
	return 1 - reduce(lambda x, y: x*y, map(lambda i: dicti[('r'+sgRNA[i]+':d'+target[i], i+1)] if sgRNA[i] != target[i] else 1, [j for j in range(0, 20)]))
	#multiple all of this in one frase. didn't did it yet. it is calleed reduce

	


#######################################################################
#
#       CRISTA: for run from any other location on the cluster
#
#######################################################################


def CRISTA(sgseq, extended29_genomic_seq, CFD_dict = None):
	"""
	predict_cleavage_score
	:param sgseq: 20 nucleotides sgRNA sequence
	:param extended29_genomic_seq: 3nt upstream + 23-nt target site + 3-nt downstream
	:return: cleavage score by CRISTA
	"""
	rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, \
	strand, start_position, end_position, include_genomic_features, w_flanking = \
		check_input_vars(rna_seq=sgseq, extended100_genomic_seq=extended29_genomic_seq,
						 genome_database=None, cell_type=None,
						 chromosome=None, strand=None,
						 start_position=None)

	features = get_features(rna_seq, dna_seq, extended100_genomic_seq, genome_database, cell_type, chromosome, strand,
							start_position, end_position, include_genomic_features, w_flanking=w_flanking)
	features_mat = np.asmatrix(features)
	score_df = predict_on_df(features_mat, include_genomic_features, w_flanking)
	return score_df.iloc[0]["CRISTA score"]
	


def find_dist(p1, p2):
	return (sum([(p1[i] - p2[i])**2 for i in range(len(p1))]))**0.5

def find_dist_np(t1, t2, p=2):
	'''p=1: City Block distance; p=2: Euclidian distance'''
	#if Euclidian:
	#	return np.linalg.norm(p1 - p2)
	return minkowski(t1, t2, p)

def find_dist_t(t1, t2, cfd_dict = None):
	p1, p2 = pos_in_metric_cfd(t1, cfd_dict), pos_in_metric_cfd(t2, cfd_dict)
	return find_dist(p1, p2)

def test_cfd_funct():
	sg = 'CCGTACGTACGTACGTACGG'
	t = 'ACGTACGTACGTACGTACGT'
	print(cfd_funct(sg,t))


def make_pos_dict(inpath = "D:\\Lab\\Cdata\\Relevant articles\\STable 19 FractionActive_dlfc_lookup.txt"):
	'''
	the dictionary manufacter here is sutable for comparing the match when the RNA sequence is represented as it's complement
	'''
	give_compl = lambda x: 'G' if x == 'C' else 'C' if x == 'G' else 'T' if x == 'A' else 'A'
	infile = open(inpath, 'r')
	next(infile)
	dicti = dict()
	#Nucs = ['A','C','T','G']
	#for Nuc in Nucs:
	#    dicti[Nuc] = [1 in range(20)]
	for line in infile:
		line_as_array = line.split('\t')
		line_as_array[0] = 'r' +  give_compl(line_as_array[0][1]) + line_as_array[0][2:]
		type, pos, score = line_as_array[0], line_as_array[1], line_as_array[5]
		dicti[(type, int(pos))] = float(score)
	pickle.dump(dicti, open("cfd_dict.p",'wb'))


def test_pos_in_metric():
	t1 = "ACGTACGTACGTACGTACGT"
	t2 = "ACGTACCCCCGTACGTACCC"
	#t3 = "ACGTACCCCCGTACGTACCC"
	p1, p2 = pos_in_metric_cfd(t1), pos_in_metric_cfd(t2)
	dist = find_dist(p1, p2)
	print(dist)
	#p1, p2 = pos_in_metric_cfd(t1), pos_in_metric_cfd(t2)
	print(find_dist(p2, p2))
	#print(dist)

def test_compre_to_shirans():
	Nuc = ["A", "G", "C", "T"]

	seq_lst = [ ''.join(random.choice(Nuc) for _ in range(20)) for i in range(200)]
	for s1 in seq_lst:
		for s2 in seq_lst:
			if cfd_funct(s1, s2) != UPGMA.cfd_func(s1, s2):
				print("different: ", s1, s2)
				print("scores :",cfd_funct(s1, s2), UPGMA.cfd_func(s1, s2))

def test2_compre_to_shirans():
	t1 = "ACGTACGTACGTACGTACGG"
	t2 = "GCGTACGTACGTACGTACGG"

	s1, s2 = 1-cfd_funct(t1, t2), 1- UPGMA.cfd_func(t1, t2)

	print(s1, s2)


if __name__ == "__main__":
	#make_pos_dict()
	#test_pos_in_metric()
	test2_compre_to_shirans()
	#make_pos_dict()
	#test_cfd_funct()