import Candidate
import os
import pickle
import statistics as st
import seaborn as sns
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h
	

def compare_10_best(path1, path2, outfile):
	disagreements = 0
	for dir in os.listdir(path1):
		res2path = "/".join([path2, dir, "res_in_lst.p"])
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		if not os.path.isfile(res2path) or not os.path.isfile(res1path):
			continue
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		#find best 10
		#first, open pickle
		res1 = pickle.load((open(res1path, "rb")))
		res2 = pickle.load((open(res2path, "rb")))
		len_of_res = min(2, len(res1))
		for i in range(len_of_res):
			if res1[i].seq != res2[i].seq: #add chacking of the score
				print (res1[i].seq)
				print(res2[i].seq)
				print(dir)
				print(i)
				disagreements +=1
				break
		print("disagreements: ", disagreements)

def compare_first(path1, path2, outfile):
	lst_dif = list()
	for dir in os.listdir(path1):
		res2path = "/".join([path2, dir, "res_in_lst.p"])
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		if not os.path.isfile(res2path) or not os.path.isfile(res1path):
			continue
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		#find best 10
		#first, open pickle
		res1 = pickle.load((open(res1path, "rb")))
		res2 = pickle.load((open(res2path, "rb")))
		if len(res1)<1 or len(res2) < 1:
		#	print("no candidates, lenres1:", len(res1), "lenres2: ",len(res2))
		#	print(dir)
			continue
		dif = (res1[0].cut_expectation - res2[0].cut_expectation)/res1[0].cut_expectation
		lst_dif.append(dif)
		#lst_dif.append((res1[0].cut_expectation - res2[0].cut_expectation)/res1[0].cut_expectation)
		if dif > 0.4:
			print('diff above 0.4: ', dif ,res1path)
	avg = st.mean(lst_dif)
	std = st.stdev(lst_dif)
	print("CI = ", mean_confidence_interval(lst_dif))
	print("CI_bayes = ",scipy.stats.bayes_mvs(lst_dif, alpha=0.95))
	print(path2)
	print("avg= ", avg)
	print("std= ", std)
	return lst_dif

		#len_of_res = min(2, len(res1))
def get_numbers(path):
	lst = list()
	for dir in os.listdir(path1):
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		if not os.path.isfile(res1path):
			continue
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		#find best 10
		#first, open pickle
		res1 = pickle.load((open(res1path, "rb")))
		if len(res1)<1:
			print("no candidates, lenres1:", len(res1))
			print(dir)
			continue
		lst.append(res1[0].cut_expectation)# /res1[0].cut_expectation)
		avg = st.mean(lst)
	print(path)
	std = st.stdev(lst)

	print("avg= ", avg)
	print("std= ", std)
	return lst
	
def get_numbers_80(path1, path2):
	lst_dif = list()
	for dir in os.listdir(path1):
		res2path = "/".join([path2, dir, "res_in_lst.p"])
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		if not os.path.isfile(res2path) or not os.path.isfile(res1path):
			continue
		res1path = "/".join([path1, dir, "res_in_lst.p"])
		#find best 10
		#first, open pickle
		res1 = pickle.load((open(res1path, "rb")))
		res2 = pickle.load((open(res2path, "rb")))
		if len(res1)<1 or len(res2) < 1:
			print("no candidates, lenres1:", len(res1), "lenres2: ",len(res2))
			print(dir)
			continue
		lst_dif.append(res1[0].cut_expectation)# /res1[0].cut_expectation)
		avg = st.mean(lst_dif)
	std = st.stdev(lst_dif)
	print(path2)
	print("avg= ", avg)
	print("std= ", std)
	return lst_dif
	
def graph():
	avg_l = [0.11220601354177076, 0.006948650820334275]
	std_l = [0.12387319637254317, 0.038570574515914856]

	
def make_DataFrame(path):
	'''return a pandas dataframe containing the data for the violin plot'''

	#creat a DataFrame object
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/"
	path1 = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409"
	paths = [path + "base_size_equals_root", path +"base_size_equals_power_0_7", path +"base_size_equals_power_0_9", path +"base_size_equals_full_500_samples"]
	x_plt = ["0.5", "0.7", "0.9", "1; sumple size = 500"] #to change
	l = []
	for p2 in paths: #make sure it is going in the corect order
		print(p2)
		l.append(compare_first(path1, p2, None))
		#l.append(get_numbers(p2))
	#l.append(get_numbers_80("/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409", paths[3]))
	#print(l)
	#l = [[],[]] + l#why?
	#l = pickle.load(open("lst_of_lst_for_sns_size_of_s.p","rb"))
	#l = pd.DataFrame(l[0:])
	print(l)
	pickle.dump(l,open("lst_of_lst_for_sns_size_of_s_dif_from_ref.p","wb"))
	#pickle.dump(l,open("lst_of_lst_for_sns_size_of_s_absolute.p","wb"))
	sns.set_style("whitegrid")
	sns_plot = sns.boxplot(data = l)# ,col = "family size", row = "expected number of cut", margin_titles=True)
	#sns_plot.set(xlabel=('2','3','4','5','6','7','8'))
	#sns_plot.savefig("violin_plot.png")
	sns_plot.figure.savefig("size_of_s.png")

	scipy.stats.bayes_mvs(data, alpha=0.95)

if __name__ == "__main__":
	#29.12:
	#path2 = '/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409'
	
	path1 = '/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/CFD_NJ'
	path2 = '/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/CFD'
	
	print('path 1: targets dist metric; path 2: cfd dist')
	print('compre first:')
	compare_first(path1, path2, None)
	exit()
	print('compre ten:')
	#compare_10_best(path1, path2, None)
	#last time I've used it:
	#path2 = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/base_size_equals_full"
	path1 = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409"
	path2 = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/base_size_equals_power_0_9"
	#path2 = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/base_size_equals_root"
	outfile = None
	#compare_first(path1, path2, outfile)
	make_DataFrame("/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/")
	
