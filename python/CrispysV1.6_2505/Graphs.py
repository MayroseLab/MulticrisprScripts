


#import pandas as pd
import pickle
import os
import Candidate
import seaborn as sns

def make_DataFrame(path):
	'''return a pandas dataframe containing the data for the violin plot'''
	#first, crate a dict. kay: size of fam. val: array of expected num of cut
	#groups_without_output = open()
	d = {}
	#d = [list()]*11
	#print(d)
	for i in range(2,11):
	#for i in range(11):

		d[i] = list()
	for dir in os.listdir(path):
		#print(dir)
		if not os.path.isfile("/".join([path, dir, "res_in_lst.p"])):
			print("asdfsadfdsafsadfsadf ", dir)
			continue
		candidates_path = "/".join([path, dir,"res_in_lst.p"])
		genes_path = "/".join([path, dir,"genesNames.p"])
		genes_lst = pickle.load(open(genes_path, "rb"))
		candidates_lst = pickle.load(open(candidates_path, 'rb'))
		if len(candidates_lst) == 0:
			print("no candidates here:")
			print(dir)
			continue
		d[len(genes_lst)].append(candidates_lst[0].cut_expectation)
		#dicti[len(genes_lst)] = dicti[len(genes_lst)] + [candidates_lst[0].cut_expectation]
	#print(d)
	#creat a DataFrame object
	#pd.DataFrame(d)
	x_plt = [str(i) for i in range(2,11)]
	l = []
	for lst in d.values(): #make sure it is going in the corect order
		l.append(lst)
	#print(l)
	l = [[],[]] + l
	pickle.dump(l,open("lst_of_lst_for_sns.p","wb"))
	sns.set_style("whitegrid")
	sns_plot = sns.violinplot(data = l)# ,col = "family size", row = "expected number of cut", margin_titles=True)
	#sns_plot.set(xlabel=('2','3','4','5','6','7','8'))
	#sns_plot.savefig("violin_plot.png")
	sns_plot.figure.savefig("output.png")
if __name__ == "__main__":
	#path = "D:\Lab\Test2"
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6/"

	make_DataFrame(path)



