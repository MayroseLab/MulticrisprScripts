__author__ = 'Gal Hyams'
##difference form the regular MC: make genes:sg dictionary in additon to the sg:genes dictionary##
import re
import os
import sys
import bottemsUpAlgorithm
import CasSites
import WrapperVers3 as wrapper
import WrapperVers3E as wrapperE
import copy
import UPGMA
import timeit
import pickle
import Metric
import make_tree_display
import set_cover_greedy
import make_tree_display_CSV
import subgroup_res

import argparse

#os.system('module load mafft/mafft7149')
#import Covers


def CRISPys_main(fasta_file, path , alg = 'A', where_in_gene = 1, use_thr = 0,  Omega = 1, df_targets = Metric.cfd_funct, protodist_outfile = "outfile", min_length= 20, max_length = 20,start_with_G = False, internal_node_candidates = 10, PS_number = 12, PAMs = 0, off_targets = 0):
	'''in dirp will be a list files. In each there will be a sequences in FASTA format. #pylip_temps_path = "/bioseq/data/results/multicrispr/test_res ## with the ability to choose where in the gene, and to give a few introns as input"'''
	start = timeit.default_timer()
	cfd_dict = None
	#print('use_thr', use_thr)
	#if isinstance(where_in_gene, str):
	#	where_in_gene = float(where_in_gene.strip())
	#if isinstance(Omega, str):
	#	Omega = float(Omega.strip())
	#if isinstance(use_thr, str):
	#	use_thr = int(use_thr.strip())
	#if isinstance(min_length, str):
	#	min_length = int(min_length.strip())
	#if isinstance(max_length, str):
	#	max_length = int(max_length.strip())
	#if isinstance(internal_node_candidates, str):
	#	internal_node_candidates = int(internal_node_candidates.strip())
	#if isinstance(PS_number, str):
	#	PS_number = int(PS_number.strip())
	#if isinstance(PAMs, str):
	#	PAMs = int(PAMs.strip())
	#df_targets = UPGMA.shalem_score
	#df_targets = UPGMA.cfd_func
	#choosing the distance function
	if df_targets == "MITScore": 
		df_targets = UPGMA.MITScore
	if df_targets == "cfd_funct" or df_targets == Metric.cfd_funct:
		df_targets = Metric.cfd_funct
		cfd_dict = pickle.load(open("/groups/itay_mayrose/galhyams/CrispysV1.6_2505/cfd_dict.p",'rb'))
	if df_targets == "ccTop":
		df_targets = UPGMA.ccTop
	if df_targets =='CRISTA':
		df_targets = Metric.CRISTA
	#df_targets = Metric.CRISTA
	if PAMs == 0:
		PAMs = ['GG']
	elif PAMs == 1:
		PAMs = ['GG','AG']
	#df_targets = UPGMA.p_distance
	#print(df_targets)
	protodist_outfile = path + "/" + protodist_outfile
	original_range_in_gene = [0, where_in_gene]
	sg_genes_dict, genes_sg_dict = dict(), dict()
	genesNames, genesList = list(), list()
	f = open(fasta_file,'r')
	gene_name = ""
	gene_seq = ""
	lines = f.readlines()
	i = 0
	genes_exons_dict = {}  #key: gene name. value: list of exons
	print()
	while i <= len(lines):
	#stage 1: make  gene: sequence dictionary
		if i == len(lines) or lines[i][0] == '>':
			if len(gene_seq) > 0 and gene_name != "": #add the gene
				if gene_name not in genes_exons_dict:
					genes_exons_dict[gene_name] = [gene_seq]
				else:
					genes_exons_dict[gene_name] = genes_exons_dict[gene_name] + [gene_seq]
				gene_seq = ""
			if i != len(lines): # lines[i-1][0] == '>':
				gene_name = lines[i][1:].strip() #without the '>' and the '\n'
		elif lines[i] != "\n":
			gene_seq += lines[i].strip()
		i+=1
	#stage 2: find the target sites
	for gene_name in genes_exons_dict.keys():
		genes_sg_dict[gene_name] = CasSites.get_targets_sites_from_exons_lst(genes_exons_dict[gene_name], original_range_in_gene, min_length, max_length,start_with_G, df_targets == Metric.CRISTA , PAMs)

		genesNames.append(gene_name)
		genesList.append("".join(genes_exons_dict[gene_name]))

		#filling up the sg_genes_dict
		for sg in genes_sg_dict[gene_name]:
			if sg in sg_genes_dict:
				sg_genes_dict[sg] = sg_genes_dict[sg] + [gene_name]
			else:
				sg_genes_dict[sg] = [gene_name]
	if alg == 'E':
		res = wrapperE.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protodist_outfile, path, df_targets, internal_node_candidates, cfd_dict, PS_number)
	else:
		res = wrapper.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protodist_outfile, path, df_targets, cfd_dict, PS_number) #thies line have been change to be sutable for wrapper
	if use_thr > 0:
		sort_thr(res, Omega, alg == 'E')
	else:
		sort_expectation(res, alg == 'E')
	pickle.dump(genesNames, open(path + "/genesNames.p", "wb"))
	pickle.dump(genesList, open(path + '/genesList.p', 'wb'))
	if alg == 'A' and use_thr > 0:
		print('use thr: ', use_thr)
		greedy_cover = set_cover_greedy.find_set_cover(res, sg_genes_dict, Omega)
		pickle.dump(greedy_cover, open(path + '/greedy_cover.p','wb'))
		make_tree_display.tree_display(path, alg == 'E', 'greedy_set_cover')
		
	#remove_repetitions_in_targets_sites(res, alg =='E')
	removed_rep = remove_repetitions_in_targets_sites(res, alg, use_thr, Omega)
	print('len res', len(res))
	print('len removed rap',len(removed_rep))
	removed_rep = removed_rep[:min(200, len(removed_rep))]
	res = res[:min(200, len(res))]
	#if len(res)>200:
	#	res = res[:200]
	wrapper.print_res_to_csvV3(res, sg_genes_dict, genesList, genesNames, path, alg =='E')

	pickle.dump(res, open(path + "/res_in_lst.p", "wb"))
	#print('removed_rep', removed_rep)
	#pickle.dump(res, open(path + '/res_in_lst_removed_rep.p','wb'))
	pickle.dump(removed_rep, open(path + '/res_in_lst_removed_rep.p','wb'))
	pickle.dump(sg_genes_dict, open(path +"/sg_genes_dict.p", 'wb'))
	pickle.dump(genes_sg_dict, open(path +"/genes_sg_dict.p", 'wb'))

	#exit()
	stop = timeit.default_timer()
	print("time: ", stop - start)
	time_file = open("time.txt", 'w')
	time_file.write(str(stop - start))
	time_file.close()
	make_tree_display_CSV.tree_display(path ,alg == 'E')
	make_tree_display.tree_display(path,alg == 'E')
	make_tree_display.tree_display(path,alg == 'E', 'removed_repetitions')
	#print('greedy cover: ',greedy_cover)
	#print(res)
	print('all well')
	return res

def sort_expectation(candidates_DS, homology):
	def sort_subgroup(candidates_DS):
		candidates_DS.sort(key = lambda item: (item.cut_expectation, item.total_num_of_mm()), reverse=True)
	if not homology:
		sort_subgroup(candidates_DS)
	else:
		for i in range(len(candidates_DS)):
			sort_subgroup(candidates_DS[i].candidate_lst)
	
def sort_thr(candidates_DS, Omega, homology):
	'''dort the candidates DS by num of genes with cut prob> Omega and then by the probobility to cleave all of these genes'''
	def sort_subgroup(candidates_DS, Omega):
		for candidate in candidates_DS:
			num_of_genes_above_thr = 0
			cleave_all = 1
			for gene, score in candidate.genes_score_dict.items():
				if score >= Omega:
					cleave_all *= score
					num_of_genes_above_thr += 1
			candidate.cleve_all_above_thr = cleave_all
			candidate.num_of_genes_above_thr = num_of_genes_above_thr
		candidates_DS.sort(key = lambda item: (item.num_of_genes_above_thr, item.cleave_all_above_thr), reverse = True)
	if not homology:
		sort_subgroup(candidates_DS, Omega)
	else:
		for i in range(len(candidates_DS)):
			sort_subgroup(candidates_DS[i].candidate_lst, Omega)
						
def sort_thr_new(candidates_DS, Omega, homology):
	'''dort the candidates DS by num of genes with cut prob> Omega and then by the probobility to cleave all of these genes'''
	def sort_subgroup(candidates_DS, Omega):
		for i in range(candidates_DS.candidate_lst):
			num_of_genes_above_thr = 0
			cleave_all = 1
			for gene, score in candidates_DS.candidate_lst[i].genes_score_dict.items():
				if score >= Omega:
					cleave_all *= score
					num_of_genes_above_thr += 1
			candidates_DS.candidate_lst[i].cleve_all_above_thr = cleave_all
			candidates_DS.candidate_lst[i].num_of_genes_above_thr = num_of_genes_above_thr
		candidates_DS.sort(key = lambda item: (item.num_of_genes_above_thr, item.cleave_all_above_thr), reverse = True)
	if not homology:
		sort_subgroup(candidates_DS, Omega)
	else:
		for i in range(len(candidates_DS)):
			sort_subgroup(candidates_DS[i], Omega)
						
def leave_only_relevant_sgRNA(res):
	if len(res) < 1:
		return
	candidates_to_del = []
	for i in range(len(res) -1, -1,-1):
		if res[i].cut_expectation < 1:
			del res[i]
		elif i < len(res) - 1:
			for j in range(i+1,len(res)):
				if j >= len(res):
					continue
				if res[i].seq == res[j].seq: # there is no need in both of them. Who is sutable for more genes?
					if res[i].cut_expectation <= res[j].cut_expectation:
						del res[i]
					else:
						del res[j]

def find_unrepresented_genes(genesNames, res):
	'''res format: list of list. each sublist: [perm sequence, int, int, list of genes, list of list: differences for each gene'''
	unrepresented = set(genesNames)  ##all the genes in the family
	represented = set()
	for perm in res:
		for gene in perm.genes_score_dict.keys():
			represented.add(gene)
	##now the represented set contains all the represented genes
	unrepresented = unrepresented - represented
	return list(unrepresented)

def add_unrepresented(genes_not_represented, genesNames, genesList):
	if len(genes_not_represented) == 1:  #find the first
		index_in_lst = genesNames.index(genes_not_represented[0])
		gene_seq = genesList[index_in_lst]
		
def remove_repetitions_in_targets_sites(candidates_lst, alg, use_thr, Omega):
	'''haven't been tested yet
	keep only with at least 1 unique targets'''
	
	def cover_above_thr(candidate, Omega, genes_list = None):
		cover_prob, num_of_covered = 1, 0
		if not genes_list:
			genes_list = [gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] >= Omega]
		for gene in genes_list:
			if candidate.genes_score_dict.get(gene,0) >= Omega:
				cover_prob *= candidate.genes_score_dict[gene]
				num_of_covered += 1
		return cover_prob, num_of_covered, genes_list
	
	def contained1(candidate_i, candidate_j, thr):
		for key, value in candidate_i.genes_score_dict.items():
			if candidate_j.genes_score_dict.get(key, 0) < value:
				return False
		return True
		
	def contained(candidate_i, candidate_j, use_thr, Omega):
		'''dose candidate_i contained in candidate_j'''
		if use_thr: 
			cover_prob_i, num_of_covered_i, genes_list = cover_above_thr(candidate_i, Omega)
		##first, check if j has all of i targets
		for key, value in candidate_i.targets_dict.items():
			if candidate_j.genes_score_dict.get(key,0) < 0.005: #epsilon
			#if not key in candidate_j.genes_score_dict:
				return False
		if use_thr:
			cover_prob_j, num_of_covered_j, _ = cover_above_thr(candidate_j, Omega, genes_list)
			if num_of_covered_j > num_of_covered_i:
				if cover_prob_j > cover_prob_i:
					print('contained')
					return True
			return False
		else:
			return candidate_i.cut_expectation < candidate_j.cut_expectation
	
	def remove_rep_subgroup(candidates_lst, use_thr, Omega):
		res = list()
		for i in range(len(candidates_lst)):
			to_append = 1
			for j in range(len(candidates_lst)):
				if (i != j and contained(candidates_lst[i], candidates_lst[j], use_thr, Omega)):
					to_append = 0
					#print('not appended')
			if to_append:
				res.append(candidates_lst[i])
		#print(res)
		#if len(res) == 0:
		#	print('len res == 0')
		return res
	
	if alg == 'E':
		res = list()
		for i in range(len(candidates_lst)):
			subgroup_candidates_lst = remove_rep_subgroup(candidates_lst[i].candidate_lst,use_thr, Omega)
			res.append(subgroup_res.Subgroup_res(candidates_lst[i].genes_lst, subgroup_candidates_lst, candidates_lst[i].name))
		return res
	else: 
		return remove_rep_subgroup(candidates_lst, use_thr, Omega)
		
		
def remove_repetitions_in_targets_sites_inplace(res, alg):
	'''haven't been tested yet
	keep only with at least 1 unique targets'''
	def remove_rep_subgroup(res):
		to_remove = list()
		targets = set()
		for i in range(len(res)):
			remove_flag = 1
			for target_tuple in res[i].targets_dict.values():
				print(target_tuple)
				if target_tuple[0][0] not in targets: #why is it [0][0].?
					#to_remove.append(i)
					remove_flag = 0
					continue
			if remove_flag == 1:
				to_remove.append(i)
			if remove_flag == 0:
				for target_tuple in res[i].targets_dict.values():
					targets.add(target_tuple[0][0])
	#now, remove
		for index in range(len(to_remove) -1, -1,-1):
			print(index)
			print(len(res))
			del res[to_remove[index]]
			
	#def should_del(res, candidate):
	#	for i in range len(res):
	#		#is there another candidate with all if it's targets?
	#		removed_flag = 1
	#		for target_tuple in res[i].targets_dict.values():
	#			print(target_tuple)
	#			if target_tuple[0][0] not in targets: #why is it [0][0].?
	#				to_remove.append(i)
	#				removed_flag = 0
					#continue
			
			
	
	if alg == 'E':
		for i in range(len(res)):
			remove_rep_subgroup(res[i])
	else: remove_rep_subgroup(res)

def parse_arguments(parser):
	#def CRISPys_main(fasta_file, path , alg = 'A', where_in_gene = 1, use_thr = 0,  Omega = 1, df_targets = Metric.cfd_funct, protodist_outfile = "outfile", min_length= 20, max_length = 20,start_with_G = False, internal_node_candidates = 10, PS_number = 12):
	parser.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The path to the input fasta file')
	parser.add_argument('path', type=str, metavar='<path>', help='The path to the directory in which the output files will be written')
	parser.add_argument('--alg', type=str, default='A', help='Choose E for considering homology')
	parser.add_argument('--where_in_gene', type=float, default=1, help='input a number between 0 to 1 in order to ignore targets sites downstream to the corresponding gene prefix')
	parser.add_argument('--t', type=int, default= 0, help='for using sgRNA to gain maximal gaining score among all of the input genes or 1 for the maximal cleavage likelihood only among genes with score higher than the average. Default: 0.')
	parser.add_argument('--v', type=float, default= 0.45, help='the value of the threshold. A number between 0 to 1 (included). Default: 0.43')
	parser.add_argument('--s', type=str, default='cfd_funct', help='the scoring function of the targets. Optional scoring systems are: cfd_funct (default), CrisprMIT and CCtop. Additinal scoring function may be added by the user or by request.')
	parser.add_argument('--p', type=str, default='outfile', help='protDist output file name. Default: "outfile"')
	parser.add_argument('--l', type=int, default= 20, help='minimal length of the target site. Default:20')
	parser.add_argument('--m', type=int, default= 20, help = 'maximal length of the target site, Default:20')
	parser.add_argument('--g', type=int, default= 0, help='1 if the target sites are obligated to start with a G codon or 0 otherwise. Default: 0.')
	parser.add_argument('--i', type=int, default= 10, help='when choosing the consider homology option, this is the number of sgRNAs designed for each homology sub-group. Default: 10')
	parser.add_argument('--ps', type=int, default= 12, help='the maximal number of possible polymorphic sites in a target. Default: 12')
	parser.add_argument('--PAMs', type=int, default=0, help='the maximal number of possible polymorphic sites in a target. Default: 12')
	parser.add_argument('--off_targets', type=str, default='', help='the maximal number of possible polymorphic sites in a target. Default: 12')

	args = parser.parse_args()
	return args

if __name__ == "__main__":
	#CRISPys_main(*sys.argv[1:])

	parser = argparse.ArgumentParser()
	args = parse_arguments(parser)
	CRISPys_main(fasta_file = args.fasta_file,
				 path = args.path,
				 alg = args.alg,
				 where_in_gene = args.where_in_gene,
				 use_thr = args.t,
				 Omega=args.v,
				 df_targets = args.s,
				 protodist_outfile = args.p,
				 min_length=args.l,
				 max_length=args.m,
				 start_with_G = args.g,
				 internal_node_candidates=args.i,
				 PS_number = args.ps,
				 PAMs = args.PAMs,
				 off_targets = args.off_targets)
	#call_using_CasSites_server_new(*sys.argv[1:])

