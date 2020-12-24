import Candidate
import Group


def fuzzy_set_cover(res, genes_lst ,size_of_output):
	'''assume res is sorted by score
	import functools
	#current_pi = functools.reduce(lambda x, y: x*y, list(map(lambda z: 1-z, lambda candidate: candidate.genes_score_dict.keys())))
	influance_on_genes = [1 for i in range(len(genes_lst))]
	objective = 0  #minimize sigma(pi(1-score)
	chosens =[res[0]] #the chosen sgRNAs
    '''


def full_cover_V1(res, genes_lst, size_of_input):
    '''find the candidates group that maximize the objective'''
    groups_DS = []
    #make the singletones
    for i in range(len(res)):
        groups_DS.append(Group(res[i].genes_score_dict, set([i])))
    #make all of the pairs



def calculate_score(lst_of_candidates, lst, genes_lst):
    '''
    :param lst_of_candidates: lst of all the candidates. "res" of algorithm A.
    :param lst: list of candidates of the current group - each candidate represented as an index
    :param genes_lst: list of the genes of the family
    :return: the objective for this group: pi sigma 1-phi(sg_j, gene)*Xj
    '''
    score = 0
    #print(len(lst_of_candidates), lst)
    for gene in genes_lst:
        gene_score = 1
        for candidate in lst:
            if gene in lst_of_candidates[candidate].genes_score_dict:
                gene_score *= (1- lst_of_candidates[candidate].genes_score_dict[gene])
        gene_score = 1- gene_score
        score += gene_score
    return score

def test_score():
    path = "D:\Lab\\test5"


def full_cover_V0(lst_of_candidates, genes_lst, k):
    '''
    Going over all of the k-mer in O(1) space. Finds the best scored one
    :param lst_of_candidates: the res from the algorithm before the cover
    :param k: size of wanted group
    :return:
    '''
    num_of_candidates = len(lst_of_candidates)
    lst = [i for i in range(k)] #thats the first set
    not_done = True
    #best_score = len(genes_lst) #find the lowest score
    best_score = 0
    while(not_done):
        score = calculate_score(lst_of_candidates, lst, genes_lst)
        #update best score and best group
        if score > best_score:
            best_score, best_group = score, lst
        not_done = increment(lst, num_of_candidates)
    return best_score, best_group


    #for i in range(k):
    #    for j in range(i,k):
     #       for t in range(j,k):
       #         for l in range(t,k):
      #              c_t = [i,j,t,l]
        #            score = score
         #           upate_best_score = score


def upper_bound(lst, index, num_of_candidates):
    return lst[index] == num_of_candidates - 1 - (len(lst)-1 -index)

def reset(lst, index):
    '''
    reset all of the digit from index and forword.
    :param lst:
    :param index: > 0
    :return:
    '''
    #print("reset")
    #print(lst,index)
    lst[index -1] += 1
    #lst[index] = lst[index-1] + 1
    #print("after reset")
    #print(lst)
    for j in range(index, len(lst)):
        lst[j] = lst[j-1] + 1
    return True


def increment(lst, num_of_candidates):
    i = len(lst) -1
    while upper_bound(lst, i, num_of_candidates):
        if i == 0:
            return False
        i -= 1
    if i == len(lst) -1:
        lst[len(lst) -1] = lst[len(lst) -1] + 1
        return True
    return reset(lst, i+1)


def increment_old(lst, num_of_candidates):
    for i in range(len(lst), -1, -1):
        if upper_bound(lst, i, num_of_candidates):
            if i == 0:
                return False
            reset(lst, i)
        else:
            lst[len(lst)] = lst[len(lst)] + 1
    return True


##################################################################################


def find_best_group_ln(res):

	'''
	import numpy as np

	scores = np.array(list(map(lambda candidate: np.array(list(candidate.genes_score_dict.values(),res)))))
	print(scores)
	scores = np.vectorize(np.log2(scores))
	print(scores)
	func = np.vectorize(lambda score: 1-score)

	scores = func(scores)
	scores = np.sum(scores, axis=1)
	print(scores)
	'''

	#rewrite after all is good
	import math

	scores = list(map(lambda candidate: list(candidate.genes_score_dict.values()),res)) #the score for each sgRNA
	print(scores)
	scores_2 = []
	i = 0
	for sg_lst_of_values in scores: #lst_of_values: the cleaving propencity for each gene
		#print(lst_of_values)

		#sg_lst_of_values = list(map(lambda x: 1 - math.log(x) if x > 0 else 0, sg_lst_of_values)) #x is the cleaving propensity of a single gene
		sg_lst_of_values = list(map(lambda x: math.log(1-x) if x < 1 else math.log(0.0001), sg_lst_of_values)) #x is the cleaving propensity of a single gene

		scores_2.append(sum(sg_lst_of_values))
		res[i].score_2 = scores_2[i]
		i += 1
	print(scores_2)
	res.sort(key= lambda candidate: candidate.score_2)

'''
    ##or, fill M in a naive way:
    M = np.zeros(len(candidate_list), len(genes_list))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            fi = 1
            for target in genes_list[j]:
                fi = fi*(1-df(candidate_list[i].seq, target)) #df gives us the cut prob for each target
        M[i][j] = 1 - fi #the probability of the j gene to be cut by the i sgRNA. In the article, the sgRNA is j, and the gene is i.

    X = new(len(candidate_list), bool) #the candidates. if x[i] == 1, the i'th candidate is chosen
    m.constrain(X.sum()) <= res_size

    #Version 1:
    for j in range(len(t)): # for each gene
        t[j] = 1- reduce(lambda x,y: x*y,map(lambda i: 1 - M[i][j]*X[i], M.shape[0]))  #there is a problem here - it will be avaluated to a number, not stay with the variables. maybe a simpler way of using lamda experssion will work?
	'''
