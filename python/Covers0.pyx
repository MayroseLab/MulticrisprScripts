import Candidate
import Group
import pickle
import copy
import os

def calculate_score(lst_of_candidates, lst, genes_lst):
    '''
    :param lst_of_candidates: lst of all the candidates. "res" of algorithm A.
    :param lst: list of candidates of the current group - each candidate represented as an index
    :param genes_lst: list of the genes of the family
    :return: the objective for this group: pi sigma 1-phi(sg_j, gene)*Xj
    '''
    cdef int score, gene_score
    score = 0
    #print(len(lst_of_candidates), lst)
    for gene in genes_lst:
        gene_score = 1
        for candidate in lst:
            if gene in lst_of_candidates[candidate].genes_score_dict:
                gene_score *= (1- lst_of_candidates[candidate].genes_score_dict[gene])
        #gene_score = 1- gene_score
        score += gene_score
    return score

def test_score():
    path = "D:\Lab\\test5"
    genes_lst = pickle.load(open(path + "\\genesNames.p", 'rb'))
    #print(genes_lst)
    lst_of_candidates = pickle.load(open(path + "\\res_in_lst.p", 'rb'))
    r = calculate_score(lst_of_candidates, [0], genes_lst)
    print(r)

def test_fc(k):
    #cdef int k
    #k = l
    path = "D:\Lab\\test5"
    genes_lst = pickle.load(open(path + "\\genesNames (2).p", 'rb'))
    #print(genes_lst)
    lst_of_candidates = pickle.load(open(path + "\\res_in_lst (2).p", 'rb'))
    r = full_cover_V0(lst_of_candidates, genes_lst, k)
    print(r)



def full_cover_V0(lst_of_candidates, genes_lst, k):
    '''
    Going over all of the k-mer in O(1) space. Finds the best scored one
    :param lst_of_candidates: the res from the algorithm before the cover
    :param k: size of wanted group
    :return:
    '''
    cdef int num_of_candidates, best_score, score, not_done
    #cdef int lst[k]
    num_of_candidates = len(lst_of_candidates)
    lst = [i for i in range(k)] #thats the first set
    not_done = True
    #best_score = len(genes_lst) #find the lowest score
    best_score = len(genes_lst)
    best_group = []
    #best_score = 0
    while(not_done):
        score = calculate_score(lst_of_candidates, lst, genes_lst)
        #update best score and best group
        if score < best_score:
            best_score, best_group = copy.deepcopy(score), copy.deepcopy(lst)
        not_done = increment(lst, num_of_candidates)
    return best_score, best_group



def upper_bound(lst, index, num_of_candidates):
    return lst[index] == num_of_candidates - 1 - (len(lst)-1 -index)

def reset(lst, index):
    '''
    reset all of the digit from index and forword.
    :param lst:
    :param index: > 0
    :return:
    '''
    cdef int j
    lst[index -1] += 1
    for j in range(index, len(lst)):
        lst[j] = lst[j-1] + 1
    return 1 #1 is True


def increment(lst, num_of_candidates):
    cdef int i
    i = len(lst) -1
    while upper_bound(lst, i, num_of_candidates):
        if i == 0:
            return 0
        i -= 1
    if i == len(lst) -1:
        lst[len(lst) -1] = lst[len(lst) -1] + 1
        return 1
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
## gready approximation ##
##################################################################################

def gready_cover(lst_of_candidates, genes_lst, k):
    '''aproximate minimal sigma pi(1-cleaving_prob(candidate, gene)
    lst_of_candidates: sorted
    '''
    group = []
    best_score = len(genes_lst)
    for n in range(k):
        #find the that minimize the score, no repetitions.
        for candidate in range(len(lst_of_candidates)):
            if (not candidate in group):
                #break
                #continue
                #print(group)
                temp_score = calculate_score(lst_of_candidates, group + [candidate], genes_lst)
                if temp_score < best_score: #update
                    best_score = copy.copy(temp_score) #in the end, the best score will be updated here
                    temp_best_candidate = copy.copy(candidate)
        if temp_best_candidate not in group:
            group.append(temp_best_candidate)
            best_score = best_score
        if best_score == 0.0:
            break
    return best_score, group


def test_gready_cover(k):
    path = "D:\Lab\\test5"
    genes_lst = pickle.load(open(path + "\\genesNames (2).p", 'rb'))
    print(genes_lst)
    lst_of_candidates = pickle.load(open(path + "\\res_in_lst (2).p", 'rb'))
    #print(lst_of_candidates)
    r = gready_cover(lst_of_candidates, genes_lst, k)
    print(r)

def run_gready_for_all():
    path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_5_11_PS"
    for dir in os.listdir(path):
        genes_path = "/".join([path, dir, "genesNames.p"])
        candidates_path = "/".join([path, dir, "res_in_lst.p"])
        f1_path =  "/".join([path, dir, "approximated_cover.txt"])
        f2_path = "/".join([path, dir, "exact_cover.txt"])

        genes_lst = pickle.load(open(genes_path, "rb"))
        candidates_lst = pickle.load(open(candidates_path, 'rb'))
        ##aproximated cover##
        f_approximated = open(f1_path, 'w')
        f_approximated.write("cover_size; cover score; cover\n")

        for i in range(2,6):

            approximated_score, approximated_cover = gready_cover(candidates_lst, genes_lst, i)
            f_approximated.write(str(i) +"; "+ str(approximated_score) + "; " + str(approximated_cover)+"\n")
        f_approximated.close()


#run_gready_for_all()



