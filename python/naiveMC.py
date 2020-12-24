__author__ = 'GH'

#the naive algorithm - for a given family, with it's group of sgRNA, find the single sgRNA that cut the maximum amount of genes
# checked in 5/7. works good, except: sg_cut_prob = 1/(distance_perm_sg+1): need to find a better way to calculat this, using shiran toll eventually

import UPGMA
import Metric
import useMUSCLE
import bisect
import copy
import Candidate
import bottemsUpAlgorithm as BU
import random

def find_Uno_sgRNA(genes_sg_dict, Omega, df, node, cfd_dict = None, PS_number = 12):
	''' not Uno any more....
	genes_sg_dict: keys are genes names, values are lists of sgRNA sequences sutable to this gene. clased = a list of the candidates DS of the children of the node '''
	##stage one: make a list of all the sgRNAs##
	lst_of_targets = list()
	for key, val in genes_sg_dict.items():
		lst_of_targets += val
	##stage two: find the sutable sgRNA:
	for_single_gene = False
	if df == Metric.CRISTA:
		initial_seq = lst_of_targets[0][:3:-6]
	else:
		initial_seq = lst_of_targets[0]
	return return_candidates(lst_of_targets, initial_seq, genes_sg_dict, Omega, df, node, for_single_gene, cfd_dict)

def return_candidates(list_of_targets, initial_seq, genes_sg_dict, Omega, df, node, for_single_gene = False, cfd_dict = None, PS_number = 12):
	dict_of_different_places = wheres_the_differences_linear(list_of_targets, df == Metric.CRISTA) ##node_targets_DS is a python array. where_the_differences.
	node.polymorphic_sites = dict_of_different_places
	#list_of_different_places = list(node.polymorphic_sites)
	if len(dict_of_different_places) > 12 :
		return None
	list_of_different_places = list(dict_of_different_places.items())
	list_of_different_places.sort(key=lambda item: item[0])
	##going over all the permutations
	list_of_perms_sequs = all_perms(initial_seq, None, list_of_different_places)
	perm_grades = []  #a list of tuples: (candidate_str,fraction_of_cut, cut_expectation, genes_list)
	#find for all permutation which genes it cut
	widest_perm_prob = 0  ## for the stopping condition. the permutation with the highst propobility to cut all the sgRNA, without considering Omega.
	lowest_of_widest_perm = 0  ## the probability the widest candidate_str will cut in the lowest cut probability sgRNA for this widest candidate_str.
	for candidate_str in list_of_perms_sequs:
		targets_dict = {} # a list of tuples: (gene name, list of target of this gene that might be cut by the candidate_str)
		wide_perm_prob = 1
		lowest_of_wide_perm = 1
		genes_covering = []  #a list of tuples: (gene name, probability to be cut).
		for gene, targets_lst_of_gene in genes_sg_dict.items(): ##find out if this gene i couched by the sgRNA seq
			prob_gene_will_not_cut = 1  ##eazier to calculate
			list_of_targets = []  ##for later knowing where the candidate_str might cut in each gene
			num_of_cuts_per_gene = 0 #in use only in the single gene version
			for target in targets_lst_of_gene:  ##targets_lst_of_gene: list of the target of the gene
				if df == Metric.CRISTA:
					#distance_candidate_target = df(candidate_str[3:-6], target)
					distance_candidate_target = df(candidate_str[3:-6], target)

				else:
					distance_candidate_target = df(candidate_str, target, cfd_dict)
				##the old one: candidate_cut_prob = 1/(distance_candidate_target+1) ##assuming distance of 1 is 100% cut. a lot of heuristics in this line: to be changed
				if distance_candidate_target == 1: #this line was changed resenetly.. it used to be "if distance_candidate_target == 0:"
					continue
				candidate_cut_prob = 1.0 - distance_candidate_target ##the distance is between 0 to 1. 0 is usually a perfect match, 1 is far
				#real sg_prob : left for later
				#change this line
				#if candidate_cut_prob > Omega /5: #want to add only the sagnificants target to this list. might need to change the thr. should be a thr that is not higher than the on in the stopping condition
				#if not BU.stopping_condition(candidate_cut_prob, Omega):
				sg_site_differents = two_sequs_differeces(candidate_str, target)
				#perm_grades[k].targets_dict[key][j] = [perm_grades[k].targets_dict[key][j], sg_site_differents]
				list_of_targets.append([target, sg_site_differents])
				prob_gene_will_not_cut = prob_gene_will_not_cut * (1- candidate_cut_prob)  #lowering the not cut prob in each sgRNA
				num_of_cuts_per_gene += candidate_cut_prob
				#for each permutation, find the probability to cut in the gene with the lowest probability to be cut
				#if lowest_of_wide_perm > candidate_cut_prob:
				#    lowest_of_wide_perm = candidate_cut_prob

			prob_gene_cut = 1 - prob_gene_will_not_cut
			if len(list_of_targets) > 0:
				targets_dict[gene] = list_of_targets  #targets of this gene to be cleaved by the current candidate


			if (for_single_gene):
				genes_covering.append((gene, num_of_cuts_per_gene))
			else:
				genes_covering.append((gene, prob_gene_cut))
			#wide_perm_prob *= prob_gene_cut  ##for the lowest of widest


		#num_of_cut = 0
		#wont_cut_prob = 1  #the probability the permutationed sequence will not cut all of the genes, that the probability each of them will be cut is greater then Omega
		cut_expection = 0.0  ##the probability the permutationed sequence will cut all of the genes, that the probability each of them will be cut is greater then Omega
		#genes_list = []  # a list of genes considered cut by this sequence
		genes_score_dict = {}  # a dict of genes: genes considered cut by this sequence, and cut prob
		for tuple in genes_covering:  #tuple : (gene name, probability to be cut)
			#if tuple[1] >= Omega:
			#num_of_cut += 1
			cut_expection += tuple[1]  ## the prob to cut all the genes
			genes_score_dict[tuple[0]] = tuple[1]
		#cut_expectation = 1 - wont_cut_prob
		#fraction_of_cut = num_of_cut/len(genes_sg_dict)  #len(genes_sg_dict) == num of genes

		##updating the targets dict##
		#for key in perm_grades[k].targets_dict.keys():  #for each gene in the the targets_dict - a list of lists. each sub list: [gene, list_of_targets]
		#    for j in range(len(perm_grades[k].targets_dict[key])): #for any match site of this gene from this specipic target. usually there is only one gene
		 #       sg_site_differents = two_sequs_differeces(perm_grades[k].seq, perm_grades[k].targets_dict[key][j])
		  #      perm_grades[k].targets_dict[key][j] = [perm_grades[k].targets_dict[key][j], sg_site_differents]

		if cut_expection >= 1: #is this condition necessary?
			current_candidate = Candidate.Candidate(candidate_str, cut_expection, genes_score_dict, targets_dict)
			perm_grades.append(current_candidate)


	del list_of_perms_sequs
	#best_perms_DS = find_max(perm_grades, genes_sg_dict)  #after the set cover
	##for finding where are the differences between the target and the DNA:

	#for k in range(len(perm_grades)): # a "best candidate_str" looks like this: (max_seq, max_fraction, max_cut_prob, genes_list, targets_dict])
	#    for key in perm_grades[k].targets_dict.keys():  #for each gene in the the targets_dict - a list of lists. each sub list: [gene, list_of_targets]
	 #       for j in range(len(perm_grades[k].targets_dict[key])): #for any match site of this gene from this specipic target. usually there is only one gene
	  #          sg_site_differents = two_sequs_differeces(perm_grades[k].seq, perm_grades[k].targets_dict[key][j])
	   #         perm_grades[k].targets_dict[key][j] = [perm_grades[k].targets_dict[key][j], sg_site_differents]
	#return perm_grades, lowest_of_widest_perm, set()
	#print(perm_grades)
	return perm_grades




def find_Uno_sgRNA_bottems_up_not_num_of_PS_stoppes(genes_sg_dict, Omega, df, node):
	''' not Uno any more....
	genes_sg_dict: keys are genes names, values are lists of sgRNA sequences sutable to this gene. clased = a list of the candidates DS of the children of the node '''
	##stage one: make a list of all the sgRNAs##
	print("node ",node)
	print(node, type(node))
	for clade in node.clades:
		print("clade: ", type(clade))
	list_of_sg = []
	for key, val in genes_sg_dict.items():
		list_of_sg += val
	##stage two: find the sutable sgRNA:
	#candidates_DS, list_of_diffrences = marge_children_candidates_DS(genes_sg_dict, Omega, df, node)
   # print("marge_children_candidates_DS")
	#print(node.clades)
	if (len(node.clades) < 1 or node.clades[0] == None) or (len(node.clades) > 1 and node.clades[1] == None):  #after making sure no bags in case only one children is present, might change this condition to "and"
		temp_candidates_DS_lst, temp_lowest_cut_site_prob = return_candidates(list_of_sg, (list(node.candidates_DS.keys()))[0], genes_sg_dict, Omega, df, node)
		temp_candidates_DS = make_candidates_dict(temp_candidates_DS_lst)
		node.candidates_DS, node.lowest_cut_site_prob = temp_candidates_DS, temp_lowest_cut_site_prob
	else:
		temp_candidates_DS, temp_lowest_cut_site_prob = marge_children_candidates_DS(genes_sg_dict, Omega, df, node)
	if not temp_candidates_DS:
		return temp_candidates_DS, temp_lowest_cut_site_prob
	#now, find the new sites which are polymorphic in the current node, but not in sub nodes:
	#for this, we will find all the polymorphic sites, of each sub node, and the polimorphic sites of the current node, and compare them.
	#this will be done using dynamic programing, in order to avoind going over all the candidates.
	#this line was written afret I came back from Italy: Did the idea was to save in each node the "list of different places", and doing BU over this as well?
   #update_node_polymorphic_sites(node)
	lowest_of_widest = find_lowest_of_widest(node)
	#find_remaining_candidates(candidate_DS, list_of_diffrences)
	#evaluate known candidates
	#print(node.candidates_DS)
	return node.candidates_DS, lowest_of_widest
	#return return_candidates(list_of_sg, list_of_sg[0], genes_sg_dict, Omega, df)

def make_candidates_dict(candidates_list):
	'''takes a list of candidate, and return a dictionary of candidates, with the seq as the key'''
	res = dict()
	for candidate in candidates_list:
		res[candidate.seq] = candidate
	return res

def marge_children_candidates_DS(genes_sg_dict, Omega, df, node):
	'''
	:param genes_sg_dict:
	:param Omega: the thr of the cut probability
	:param df: distance function
	:param clades_candidates_DS: a data structure contains the candidates of the clade. cl
	:return:
	'''
	node.set_candidates_DS(node.candidates_DS)
	#node.set_candidates_DS()
	candidates_DS = node.candidates_DS
	#print("first time print candidate ds, clades[1]; ", node.clades[1].candidates_DS)
	children_candidates_DS = list(map(lambda x:x.candidates_DS, node.clades))

	#print("children C DS:\n",children_candidates_DS)
	#print(children_candidates_DS)
	for i in range(len(children_candidates_DS)): #len(clades_candidates_DS) == 2. "for each child node"

		if (children_candidates_DS[i] == None or children_candidates_DS[(i+1)%2] == None): #for now, there is no need in this
			continue
		#if children_candidates_DS[(i+1)%2] == None:
		#    other_clade_polymorphic_sites = dict()
		#else:
		#other_clade_polymorphic_sites = wheres_the_differences(list(map(lambda x:str(x), children_candidates_DS[(i+1)%2].keys()))) ##node_targets_DS is a python array. where_the_differences.
		other_clade_polymorphic_sites = node.clades[(i+1)%2].polymorphic_sites
		print(other_clade_polymorphic_sites)
			# print("other clade PS print1: ", other_clade_polymorphic_sites)
#            if len(other_clade_polymorphic_sites) > 3:
 #               print("other calde PS >5")
				#return(dict(), 0)
  #              return(None, None)

		for candidate_str in children_candidates_DS[i].keys():  # the key is this dictionary, is the sequnce of the candidate sgRNA
			#other_clade = children_candidates_DS[(i+1)%2]
			#other_clade_polymorphic_sites = list(wheres_the_differences_BU(other_clade.).items()) ##node_targets_DS is a python array. where_the_differences.
			#other_clade_polymorphic_sites = list(wheres_the_differences(children_candidates_DS[(i+1)%2].keys()).items()) ##node_targets_DS is a python array. where_the_differences.
			#print("0", list(map(lambda x:str(x), wheres_the_differences(children_candidates_DS[(i+1)%2].keys()))))
			#print("1", wheres_the_differences(children_candidates_DS[(i+1)%2].keys()))
			#other_clade_polymorphic_sites = list(wheres_the_differences(list(map(lambda x:str(x), children_candidates_DS[(i+1)%2].keys())))) ##node_targets_DS is a python array. where_the_differences.

			temp_candidates_DS = list() #for not changing the size of the in itereted dictionary druing the iteration
			for gene in children_candidates_DS[i][candidate_str].genes_score_dict.keys():
				#print("gene", gene)
				if children_candidates_DS[i][candidate_str].genes_score_dict[gene] < Omega: #shuoldn't add the candidate to the structure, when considering the current gene. this is the BNB (branch and bound)
					break
				#if we got to here, at least 1 gene is expected to be cleaved with high probability
					#first, cheack if there is redundency of this candidate in the other clade
				current_candidate = copy.deepcopy(children_candidates_DS[i][candidate_str]) # after making sure no bugs, to change to not a deep copy.
				if i == 0 and current_candidate in children_candidates_DS[(i+1)%2].keys():
					# a redundency. update and add
					#update the fraction of cut
					other_genes_score_dict = children_candidates_DS[(i+1)%2][candidate_str].genes_score_dict
					other_targets_dict = children_candidates_DS[(i+1)%2][candidate_str].match_sites_dict
					for gene_name in other_genes_score_dict.keys():
						current_candidate.add_known_site(other_genes_score_dict, other_targets_dict[gene_name], gene_name)
						#clades_candidate[i][candidate].add_known_site(self, other_genes_score_dict, other_targets_dict[key], key)
					#append_candidate_to_candidates_DS(current_candidate, candidtes_DS)
				else:  #find if it has any influance on the targets of the other node
					for other_candidate_str in children_candidates_DS[(i+1)%2].keys():
						other_genes_score_dict = children_candidates_DS[(i+1)%2][other_candidate_str].genes_score_dict
						other_targets_dict = children_candidates_DS[(i+1)%2][other_candidate_str].targets_dict
						for gene_name in other_genes_score_dict.keys():
							#children_candidates_DS[i][current_candidate.seq].add_known_site(other_genes_score_dict, other_targets_dict[gene_name], gene_name)
#                            current_candidate.add_known_site(other_genes_score_dict, other_targets_dict[gene_name], gene_name)
							current_candidate.add_known_sites_of_gene(other_targets_dict, other_genes_score_dict[gene_name], gene_name)


				#current_candidate.recalculate_cut_prob(Omega,len(genes_sg_dict)) #haven't cheacked it if this line is in the right place
				#print("append to temp C DS")
				temp_candidates_DS.append(current_candidate)

			for current_candidate in temp_candidates_DS:
				#print("call for append_candidate_to_candidates_DS")
				append_candidate_to_candidates_DS(current_candidate, node)
				#now, find the new candidate by apllying the polimorphoc sites from the other node on the current_candidate
			#is the indentation correct? does one more \t is needed?
			#print("other calde PS: ", other_clade_polymorphic_sites)
			if not isinstance(other_clade_polymorphic_sites, list):
				other_clade_polymorphic_sites = list(other_clade_polymorphic_sites.items())
			for polymorphic_site in other_clade_polymorphic_sites:
				for j in range(len(polymorphic_site[1])):
					candidets_to_append = []
					for current_candidate in candidates_DS.values():#newxt to do: swhich the order of those loops
						if current_candidate.seq[polymorphic_site[0]] != polymorphic_site[1][j]: #
							current_seq = current_candidate.seq[:polymorphic_site[0]] + polymorphic_site[1][j] + current_candidate.seq[polymorphic_site[0]+1 :]
							#check relation with the target of the sub tree
							new_current_candidate = create_a_new_candidate_and_fill_fields(current_seq, genes_sg_dict, df, Omega)
							#if not (new_current_candidate):
							#    continue
							candidets_to_append.append(new_current_candidate)
					for candidate in candidets_to_append:
						append_candidate_to_candidates_DS(candidate, node)
							#if (new_current_candidate):
							#    candidates_DS[current_seq] = new_current_candidate
	node.set_polymorphic_sites(updated_node_polymorphic_sites(node))
	return node.candidates_DS, node.lowest_cut_site_prob

def create_a_new_candidate_and_fill_fields(current_seq, genes_sg_dict, df, Omega):
	#next stage of work: add lowest_cut_site_prob
	genes_score_dict = {}
	targets_dict = {}
	number_of_node_genes = len(genes_sg_dict)
	for gene_name, genes_targets_list in genes_sg_dict.items():
		prob_gene_will_not_cut = 1  ##eazier to calculate
		list_of_targets = []  ##for later knowing where the perm might cut in each gene
		for target in genes_targets_list:  ##sg_lst_of_gene: list of the sg of the gene
			distance_candidate_target = df(current_seq, target)
				##the old one: sg_cut_prob = 1/(distance_perm_sg+1) ##assuming distance of 1 is 100% cut. a lot of heuristics in this line: to be changed
			candidate_target_cut_prob = 1 - distance_candidate_target ##the distance is between 0 to 1. 0 is usually a perfect match, 1 is far
				#real sg_prob : left for later
			# to change this line so it will be compatible with the Bottem up stopping condition
			#if candidate_target_cut_prob > Omega /5: #want to add only the sagnificants sg to this list. might need to change the thr. should be a thr that is not higher than the on in the stopping condition
			if not BU.stopping_condition(candidate_target_cut_prob, Omega):
				list_of_targets.append(target)
				prob_gene_will_not_cut = prob_gene_will_not_cut * (1- candidate_target_cut_prob)  #lowering the not cut prob in each sgRNA
				#for each permutation, find the probability to cut in the gene with the lowest probability to be cut
		prob_gene_cut = 1 - prob_gene_will_not_cut
		if prob_gene_cut > Omega:
			genes_score_dict[gene_name] = prob_gene_cut
		fraction_of_cut = len(genes_score_dict)/number_of_node_genes
		#make the match_site_dict
		if not len(list_of_targets) == 0:
			match_sites_dict_value = []
			for target_site in list_of_targets:
				match_sites_dict_value.append([target_site, two_sequs_differeces(current_seq, target_site)])
			targets_dict[gene_name] = match_sites_dict_value
	res = Candidate.Candidate(current_seq, fraction_of_cut, prob_gene_cut, genes_score_dict, targets_dict)
	return res

def append_candidate_to_candidates_DS(current_candidate, node):
	#these two lines are not true

	if not current_candidate:
		return
	if current_candidate.seq in node.candidates_DS:
		#update the exsisting one: add target sites and recalculate the scores
		node.candidates_DS[current_candidate.seq].add_known_sites(current_candidate, node)  # add update of lowest_of_widest to that method. the node is given as well, for knowing the number of genes in this node.
		#update fraction of cut and node's widest and lowest of widest
		#print("whats goin bro")
	else:
		#print("current candidate seq not in node.candidates_DS")
		#put a new Candidate in the node candidates_DS
		node.candidates_DS[current_candidate.seq] = current_candidate
		if node.lowest_cut_site_prob > current_candidate.lowest_cut_site_prob:
			node.lowest_cut_site_prob = current_candidate.lowest_cut_site_prob

	#if node.widest < current_candidate.widest:
	 #   node.lowest_cut_site_prob

def updated_node_polymorphic_sites(node):
	'''might have bugs'''
	#base condition
	for i in range(len(node.clades)):
		if node.clades[i].is_terminal: # can be done better
			return wheres_the_differences(node.node_targets_DS)
	child_a_representor, child_b_representor = random.choice(list(node.clades[0].candidates_DS.keys())), random.choice(list(node.clades[1].candidates_DS.keys()))
	#base condition:
	#i don't think this case can happen. the condition above should cover this case.
	if not child_a_representor:
		if not child_b_representor:
			return wheres_the_differences(node.node_targets_DS)
		else: #node b is a leave: contains only one candidate
			return copy.deepcopy(node.clades[1].polymorphic_sites)
	if child_a_representor and not child_b_representor:
		return copy.deepcopy(node.clades[0].polymorphic_sites) #we might don't want a deep copy here
	child_a_ps = node.clades[0].polymorphic_sites
	child_b_ps = node.clades[1].polymorphic_sites
	res = merge_dicts(child_a_ps, child_b_ps)
	#now, find the new polymorfic sites by finding the diferences between representors from each of the children
	if len(child_a_representor) > len (child_b_representor):
		longer = child_a_representor
		shorter = child_b_representor
	else:
		longer = child_b_representor
		shorter = child_a_representor
	for i in range(len(shorter)):
		if i not in res and shorter[i] != longer[i]:
			res[i] = [child_a_representor[i], child_b_representor[i]]
	for i in range(1,len(longer) - len(shorter)):
		res[len(longer) - i] = longer[len(longer) -i]
	return res

def merge_dicts(d1, d2):
	'''
	:param d1: dictionary. value is a list of items
	:param d2: the same
	:return: merge to dictionaries when the values in the dictionaries are lists
	'''
	res = dict()
	for key, value in d1.items():
		if key in d2:
			new_value = list(set(d1[key] + d2[key]))
			res[key] = new_value
		else:
			res[key] = d1[key]
	for key, value in d2.items():
		if key not in d1:
			res[key] = d2[key]
	return res


def wheres_the_differences_BU(leave_DS, known_polymorphic_sites_a, known_polymorphic_sites_b, node_a_representor, node_b_representor):
	'''
	the current challange: how to make it robust to different lengths of the targets sites?
	leave_DS: a data structure containing the leaves
	known_polymorphic_sites:
	'''
	#base condition:
	if not node_a_representor:
		if not node_b_representor:
			return wheres_the_differences(leave_DS) #this is a dictionary. key: place of diference. value: what is the diference.
		else: #node b is a leave: contains only one candidate
			return copy.deepcopy(known_polymorphic_sites_b)
	if node_a_representor and not node_b_representor:
		return copy.deepcopy(known_polymorphic_sites_a)
	else:
		res = copy.deepcopy(known_polymorphic_sites_a)
		for key, value in known_polymorphic_sites_b.items():
			if key not in res:
				res[key] = value
			else:
				res[key] = res[key] + value
		if len(node_a_representor) > len (node_b_representor):
			longer = node_a_representor
			shorter = node_b_representor
		else:
			longer = node_b_representor
			shorter = node_b_representor
		for i in range(len(shorter)):
			if i not in res and shorter[i] != longer[i]:
				res[i] = [node_a_representor[i], node_b_representor[i]]
		for i in range(1,len(longer) - len(shorter)):
			res[len(longer) - i] = longer[len(longer) -i]
				#should
		return res  #if it is a

	differences = {}  #key: place where there is a difference. value: letter apeared in this place
	for i in range(len(leave_DS)): ##node_targets_DS is a python array
		for j in range(i, len(leave_DS)):
			current_differences = two_sequs_differeces(leave_DS[i], leave_DS[j])
			for t in current_differences:
				#if t in differences:
				#   differences[t] = current_differences[t]
				#else:
				if t in differences:
					differences[t] = differences[t] + current_differences[t]
				else:
					differences[t] = current_differences[t]
	return differences

def find_lowest_of_widest(node):
	return node.lowest_cut_site_prob

def find_remaining_candidates(candidates_DS, list_of_diffrences):
	#first, find the new polimorphic sites: the
	dict_of_different_places = wheres_the_differences(list_of_sg) ##node_targets_DS is a python array. where_the_differences.
	#if len(dict_of_different_places) > 2: #change to 11
	#    return None, None
	list_of_different_places = list(dict_of_different_places.items())
	list_of_different_places.sort(key=lambda item: item[0])

def find_max(tup_lst, genes_sg_dict):
	'''this algorithm is wrong, but not in useage'''
	#initial list of genes to cover
	genes = list(genes_sg_dict.keys())
	tup_lst.sort(reverse= True, key = lambda tup: (tup[1], tup[2]))  #sort by fraction of cut genes, and then by score of cutting.
	# reminder: a "best perm" looks like this: (max_seq, max_fraction, max_cut_prob, genes_list, match_sites_dict)
	# now, add to res each relevant perm, by order.
	res = []
	for tup in tup_lst:
		for gene_name in tup[3]:
			if gene_name in genes:
				res.append(tup)
				for gene_cover_by_tup in tup[3]:
					if gene_cover_by_tup in genes:
						genes.remove(gene_cover_by_tup)
				break
	return res

def find_max_for_a_single_sg(tup_lst):
	'''input is a lst of tuples: (perm, fraction genes being cut among all the genes, probability to cut all the genes in genes list, genes_list)'''
	max_fraction = 0
	for tuple in tup_lst:
		if tuple[1] > max_fraction:
			max_fraction = tuple[1]
	##removing the unimportant permutation sequence
	i=0
	while i < (len(tup_lst)):
		if tup_lst[i][1] < max_fraction:
			del tup_lst[i]
		else:
			i+=1
	##find the best permutation sequence among all the permutation sequences that cut the maximum amount of genes
	max_cut_prob = 0  #the prob to cut all of the genes been cut: prob(cut1)*prob(cut2)*...*prob(cut_n)
	max_seq = ""
	for tuple in tup_lst:
		if tuple[2] > max_cut_prob:
			#max_fraction = tuple[1]
			max_cut_prob = tuple[2]
			res = tuple
			#max_seq = tuple[0]
			#genes_list = tuple[3]
	#return (max_seq, max_fraction, max_cut_prob, genes_list)  #number of cut genes == len(genes_list)
	return res

def find_mighty_max(tup_lst):
	'''
	just like find max, but sorting the list, and astemeting the chosen sg quality. if it is not good enoght, going to the next option.
	'''
	#find the list length
	max_fraction = 0
	fractions_avilable = []
	sublists = []  # a list of lists. each sublist contains all the tuples representing sequences with the same fraction
	for tuple in tup_lst:
		if tuple[1] > max_fraction:
			max_fraction = tuple[1]
		insertAt = bisect.bisect_left(fractions_avilable, tuple[1])
		if sublists[insertAt] != tuple[1]:
			fractions_avilable.insert(insertAt, tuple[1])
			sublists.insert(insertAt, [])
		sublists[insertAt].append(tuple)
	while sublists:
		current_sublist_index = len(sublists) -1
		sublists[current_sublist_index].sort(lambda tuple: tuple[2])
		while len(sublists[current_sublist_index]) > 0:
			candidate = sublists[current_sublist_index][-1][0]
			if test_candidate(candidate):
				return candidate
			else:
				sublists[current_sublist_index].pop()  ##remove the last item from the list
				if len(sublists[current_sublist_index]) == 0:
					sublists.pop()
	return None  #check rather is should be return None, None

def test_candidate(candidate):
	'''update: will be done at the beginning of the MC algorithm
	test the CG content, off target effect, and maybe other elements will be added'''
	num_of_CG = 0
	for i in range(len(candidate)):
		if candidate[i] == 'C' or candidate[i] == 'G':
			num_of_CG += 1
	if num_of_CG/len(candidate) < 1/5 or num_of_CG/len(candidate) > 4/5:
		return False
	#now check for off targets
	return True

def simpleTest():
	'''testing uno'''
	#genes_sg_dict = {"gene1": ["acgtc", "atttc"] , "gene2": ["atcct"], "gene3": ["actct", "accct"] }
	genes_sg_dict = {"gene1": ["acgtacgt".upper()], "gene2": ["acgtagct".upper()], "gene3": ["acgtattg".upper()], "gene4": ["atgcacgt".upper()], "gene5": ["atgcatgc".upper()]}
	print("res:", find_Uno_sgRNA(genes_sg_dict, 0.011))

def two_sequs_differeces_by_indeces(seq1,seq2):
	'''return a list of where the two sequences are different'''
	differences = {}  ##key: place of disagreement. value: the suggestions of each side
	set_of_differences = {'A', 'C','G', 'T'}
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	for i in range(1,len(seq2) - len(seq1)):
		differences[len(seq2) - i] = copy.copy(set_of_differences)
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			differences[i] = copy.copy(set_of_differences)
	return differences




###############from the old algorithm###########


def wheres_the_differences_not_naive(leave_DS):
	'''return a list of places which at least two sequences are different at'''
	checked_group = set() ##all the memmbers at this group have already been checked reagarding each other
	differences_in_checked_group = []
	for i in range(len(leave_DS)): ##node_targets_DS is a python array
		differences_from_checked_group = []
		not_to_check = set()
		##seen_cheaked = False  ##if the sequence had alrady tested against a memmber from the cheacked group
		for j in range(i, len(leave_DS)):
			if j in checked_group:
				current_differences = two_sequs_differeces(leave_DS[i], leave_DS[j], not_to_check)
				differences_from_checked_group += current_differences
				for i in current_differences:
					not_to_check.add(i)


def all_perms(initial_seq, list_of_sequs, list_of_differences):
	'''each recursive call add the next part to the sequnces. the resuls is sequnces off each of the parms
	list of differences : list of tuples: (place, set of letters)'''
	if len(list_of_differences) == 0:  ##the stopping condition
		if(list_of_sequs):
			return list_of_sequs
		elif(initial_seq):
			return [initial_seq]
		else:
			return []
	else:
		new_list_of_sequs = []
		if not (list_of_sequs):  ##initialising the list of sequences
			list_of_sequs = []
			list_of_sequs.append(initial_seq[:list_of_differences[0][0]])
		#else:
		for seq in list_of_sequs:
			for letter in list_of_differences[0][1]:
				if len(list_of_differences) > 1:
					new_list_of_sequs.append(seq + letter + initial_seq[len(seq)+1 :list_of_differences[1][0]]) #the place of the next versital letter place
				else:
					new_list_of_sequs.append(seq + letter + initial_seq[len(seq)+1 :])
		del list_of_sequs
		return all_perms(initial_seq, new_list_of_sequs, list_of_differences[1:])


def test_all_perms():
	initial_seq = "acgtacgt" #chose a leaf randomly
	leaves_DS = ["acgtacgt", "acatacgt", "acttacgt", "acgtacct"]
	list_of_differencces = list(wheres_the_differences(leaves_DS).items())
	list_of_differencces.sort(key = lambda item: item[0])
	print(list_of_differencces)
	perms = all_perms(initial_seq, None, list_of_differencces)
	print(perms)

def wheres_the_differences_specific(leave_DS):
	''' return a dict: key is a place in which at least two sequences are different at, and value is a set of each letter that was in a seq in this place '''
	differences = {}  #key: place where there is a difference. value: letter apeared in this place
	for i in range(len(leave_DS)): ##node_targets_DS is a python array
		for j in range(i, len(leave_DS)):
			current_differences = two_sequs_differeces_by_indeces(leave_DS[i], leave_DS[j])
			for t in current_differences:
				#if t in differences:
				#   differences[t] = current_differences[t]
				#else:
				if t in differences:
					differences[t] = differences[t] | current_differences[t]
				else:
					differences[t] = current_differences[t]
	return differences


def two_sequs_differeces(seq1,seq2, CRISTA = False):
	'''return a list of where the two sequences are different'''
	if CRISTA:
		seq1, seq2 = seq1[3:-6], seq2[3:-6]
	differences = {}  ##key: place of disagreement. value: the suggestions of each side
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	for i in range(1,len(seq2) - len(seq1)): #when the length of the sequences are different
		differences[len(seq2) - i] = seq2[len(seq2) -i]
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			#differences[i] = list(set([seq1[i], seq2[i]])
			differences[i] = [seq1[i], seq2[i]]
	return differences

def wheres_the_differences(leaves_DS):
	differences = {}  #key: place where there is a difference. value: letter apeared in this place
	for i in range(len(leaves_DS)): ##node_targets_DS is a python array
		for j in range(i+1, len(leaves_DS)):
			current_differences = two_sequs_differeces(leaves_DS[i], leaves_DS[j])
			for t in current_differences:
				#if t in differences:
				#   differences[t] = current_differences[t]
				#else:
				if t in differences:
					differences[t] =list(set(differences[t] + current_differences[t]))
				else:
					differences[t] = current_differences[t]
	return differences

def wheres_the_differences_linear(leaves_DS, CRISTA = False):
	differences = dict()  #key: place where there is a difference. value: letter apeared in this place
	if len(leaves_DS)<2:
		return differences
	ref = leaves_DS[0]
	for i in range(1,len(leaves_DS)): ##node_targets_DS is a python array
		current_differences = two_sequs_differeces(ref, leaves_DS[i], CRISTA)
		for t in current_differences:
				#if t in differences:
				#   differences[t] = current_differences[t]
				#else:
			if t in differences:
				differences[t] =list(set(differences[t] + current_differences[t]))
			else:
				differences[t] = current_differences[t]
	return differences


def wheres_the_differences_BU_one_child(leave_DS, known_polymorphic_sites_a, node_a_representor, node_b_representor):
	res = {}
	for i in range(len(node_a_representor)):
		if i not in known_polymorphic_sites_a:
			if node_a_representor[i] != node_b_representor[i]:
				res[i] = {}



if __name__ == "__main__":
	simpleTest()
