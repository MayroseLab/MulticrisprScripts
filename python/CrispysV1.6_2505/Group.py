import Candidate

class Group:
    def __init__(self, genes_score_dict, candidates_set):
        '''
        :param genes_score_dict: the score for each gene: pi(1-score(sg,gene).
        :param candidates_set:
        self.score: summation over all of the values in genes_score_dict
        :return:
        '''
        self.genes_score_dict = genes_score_dict or {}
        self.candidates_set = candidates_set or set()
        self.score = 0

    def add_candidate(self, candidate_number, res):
        if candidate_number in self.candidates_set:
            return
        self.candidates_set.add(candidate_number)
        for gene, score in res[candidate_number].genes_score_dict.items():
            if gene in self.genes_score_dict:
                self.genes_score_dict[gene] = self.genes_score_dict[gene]  * (1-score)
            else:
                self.genes_score_dict[gene] = score
