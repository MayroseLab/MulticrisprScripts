import Candidate
import Covers2
import os
import pickle
import numpy
import seaborn as sns
import set_cover_greedy
import pandas as pd
import numpy as np

def test_compare_greedy_regular(path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11", thr = 0.5):
	res_file = open('05_02_tests.txt', 'w')
	for dir in os.listdir(path):
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		genes_sg_path = "/".join([path, dir, "genes_sg_dict.p"])
		sg_genes_path = "/".join([path, dir, "sg_genes_dict.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		#if not best_thr_candidate:
		#	problematics_f.write(dir + ": empty candidates lst\n")
		#	continue
		#if best_amount == num_of_genes:
			#update_d(full_cover,num_of_genes,1)

		#update_d(num_of_fam, num_of_genes, 1)

		set_cover = Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path)
		sg_genes_dict, candidates_lst = pickle.load(open(sg_genes_path,'rb')), pickle.load(open(candidates_path, 'rb'))
		genes_sg_dict = pickle.load(open(genes_sg_path,'rb'))
		greedy_sc = set_cover_greedy.find_set_cover(candidates_lst, sg_genes_dict, thr, genes_sg_dict)
		if len(greedy_sc) < len(set_cover):#greedy is better:
			res_file.write(str(dir)+"\n")
			res_file.write('greedy: '+ str(greedy_sc)+'\n')
			res_file.write('regular: '+ str(set_cover)+'\n')
			for gene in genes_sg_dict.keys():
				if len(genes_sg_dict[gene]) == 0:
					res_file.write("no targets for gene " + str(gene) +"\n")
	res_file.close()

def tables_for_greedy_set_cover(path, outpath, thr, graps = False):
	if graps == False:
		sc_dict = {i:0 for i in range(2,11)}	#for the size of the set cover
	else:
		sc_dict = dict()
	num_of_fam = {i:0 for i in range(2,11)}
	for dir in os.listdir(path):
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		genes_sg_path = "/".join([path, dir, "genes_sg_dict.p"])
		sg_genes_path = "/".join([path, dir, "sg_genes_dict.p"])
		#problematics_f = open("/".join([outpath, "problematics_greedy_set_cover_31_12.txt"]),'w')
		if not os.path.isfile(candidates_path):
			#problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))	
		update_d(num_of_fam, num_of_genes, 1)
		
		#set_cover = Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path)
		sg_genes_dict, candidates_lst = pickle.load(open(sg_genes_path,'rb')), pickle.load(open(candidates_path, 'rb'))
		genes_sg_dict = pickle.load(open(genes_sg_path,'rb'))
		set_cover = set_cover_greedy.find_set_cover(candidates_lst, sg_genes_dict, thr, genes_sg_dict)
		key = num_of_genes #think that kay sepose to be num of gene family.
		if graps:
			update_d_SD(sc_dict, key, len(set_cover))
		else:
			update_d(sc_dict, key, len(set_cover))
	if graps:
		x_plt = [str(i) for i in range(2,11)]
		l = []
		for lst in sc_dict.values(): #make sure it is going in the corect order
			l.append(lst)
		#print(l)
		l = [[],[]] + l
		pickle.dump(l,open("lst_of_lst_for_sns_0_45_greedy.p","wb"))
		#sns.set_style("whitegrid")
		#sns_plot = sns.violinplot(data = l)# ,col = "family size", row = "expected number of cut", margin_titles=True)
		#sns_plot.set(xlabel=('2','3','4','5','6','7','8'))
		#sns_plot.savefig("violin_plot.png")
		#sns_plot.figure.savefig("greedy_set_cover_violin_plot_0_45.png", dpi = 600)
	pickle.dump(sc_dict,open(outpath + '/greedy_sc_dict.p', 'wb'))
	pickle.dump(num_of_fam, open(outpath + '/num_of_fum_greedy_sc_dict.p','wb'))
	with open(outpath + '/greedy_sc_res_thr_'+str(thr)+'_26_01.txt', 'w') as csvfile:
		csvfile.write("fam size;number of fam;sc\n")
		for i in range(2,11): #the sizes of the families
			if i in num_of_fam:
					#print(i)
					#print(num_of_fam[i])
					#print(sc_dict)#/num_of_fam[i]
				csvfile.write(str(i) +";"+ str(num_of_fam[i])+";"+ str(sc_dict[i]/num_of_fam[i])+'\n')
				#";"+ str(fst_dict[i]/num_of_fam[i])+ ";"+ str(resc_dict2[i]/num_of_fam[i]) + ";"+ str(resc_dict3[i]/num_of_fam[i])+";"+ str(resc_dict4[i]/num_of_fam[i])+ ";"+ str(resc_dict5[i]/num_of_fam[i]) +"\n")

def tables_of_covers(path, outpath, thr):
	'''to each gene family, finds the size of the set cover, fuzzy set cover, and restricted set cover with values of 1 to 5'''
	sc_dict = {i : 0 for i in range(2,11)}
	fst_dict = {i : 0 for i in range(2,11)}

	resc_dict2, resc_dict3, resc_dict4, resc_dict5 = {i : 0 for i in range(2,11)},  {i : 0 for i in range(2,11)},  {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}
	num_of_fam =  {i : 0 for i in range(2,11)}
	full_info_dict = dict()
	problematics_f = open(outpath + '/problematics.txt', 'w')
	for dir in os.listdir(path):
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		genes_sg_path = "/".join([path, dir, "genes_sg_dict.p"])
		sg_genes_path = "/".join([path, dir, "sg_genes_dict.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		#if not best_thr_candidate:
		#	problematics_f.write(dir + ": empty candidates lst\n")
		#	continue
		#if best_amount == num_of_genes:
			#update_d(full_cover,num_of_genes,1)
			
		update_d(num_of_fam, num_of_genes, 1)
		
		set_cover = Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path)
		print("c_path: ", candidates_path)
		f_sc = Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path)
		cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
		key = num_of_genes #think that kay sepose to be num of gene family.
		if f_sc:
			update_d(fst_dict, key, len(f_sc)) 
		if set_cover:
			update_d(sc_dict, key, len(set_cover))
		update_d(resc_dict2, key, cover2[0]) #covers_i[0] is 'best score'
		update_d(resc_dict3, key, cover3[0])
		update_d(resc_dict4, key, cover4[0])
		update_d(resc_dict5, key, cover5[0])
		full_info_dict[dir] = [resc_dict2, resc_dict3, resc_dict4, resc_dict5]
	pickle.dump(full_info_dict,open(outpath + '/full_info_covers_dict.p', 'wb'))
	with open(outpath + '/sc_res_21_01.txt', 'w') as csvfile:
		csvfile.write("fam size;number of fam;sc;fsc;resc_dict2; resc_dict3; resc_dict4; resc_dict5\n")
		for i in range(2,11): #the sizes of the families
			if i in num_of_fam and i in sc_dict and i in fst_dict and i in resc_dict2 and i in resc_dict3 and i in resc_dict4 and i in resc_dict5:
					#print(i)
					#print(num_of_fam[i])
					#print(sc_dict)#/num_of_fam[i]
				csvfile.write(str(i) +";"+ str(num_of_fam[i])+";"+ str(sc_dict[i]/num_of_fam[i])+";"+ str(fst_dict[i]/num_of_fam[i])+ ";"+ str(resc_dict2[i]/num_of_fam[i]) + ";"+ str(resc_dict3[i]/num_of_fam[i])+";"+ str(resc_dict4[i]/num_of_fam[i])+ ";"+ str(resc_dict5[i]/num_of_fam[i]) +"\n")
				#print(dir)		

def diff_versatile_thr_greedy_regular_fuzzy_covers():
	'''CFD for now
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	general_sc_dict = dict() #key:thr; val: dict of dif
	general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		sc_dict = {i : 0 for i in range(2,11)}
		fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#
			update_d(sc_dict, key, (len_set_cover - len_greedy)/len_greedy)
			update_d(fsc_dict, key, (len_f_sc - len_greedy)/len_greedy)
		#was not used to make the original file:
		fam_sizes = [i for i in range(2,11)]
		for size in fam_sizes:
			sc_dict[size] = sc_dict[size]/num_of_fam[size]
			fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		general_sc_dict[thr] = sc_dict
		general_fsc_dict[thr] = fsc_dict
	problematics_f.close()
	pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_greedy_per_thr_07_02', 'wb'))
	pickle.dump(general_sc_dict, open(path + '/sc_dict_diff_from_greedy_per_thr_07_02', 'wb'))


def diff_versatile_thr_fuzzy_vs_regular():
	'''CFD for now
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	#general_sc_dict = dict() #key:thr; val: dict of dif
	general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		#sc_dict = {i : 0 for i in range(2,11)}
		fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			#len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#update_d(sc_dict, key, (len_set_cover - len_greedy)/len_greedy)
			update_d(fsc_dict, key, (len_f_sc - len_set_cover)/len_set_cover)
		#was not used to make the original file:

		#
		fam_sizes = [i for i in range(2,10,2)]
		new_fsc_dict = dict()
		for size in fam_sizes:
			if size != 8:
				new_fsc_dict[(size, size+1)] = (fsc_dict[size] + fsc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_fsc_dict[(size, size+1, size+2)] = (fsc_dict[size] + fsc_dict[size+1] + fsc_dict[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])

		#new_fsc_dict[(10)] =
		#
		#fam_sizes = [i for i in range(2,11)]
		#for size in fam_sizes:
			#sc_dict[size] = sc_dict[size]/num_of_fam[size]
		#	fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		#general_sc_dict[thr] = sc_dict
		general_fsc_dict[thr] = new_fsc_dict
	problematics_f.close()
	pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_regular_per_thr_20_02_coupled', 'wb'))
	#pickle.dump(general_sc_dict, open(path + '/sc_dict_diff_from_greedy_per_thr_07_02', 'wb'))

def fraction_fuzzy_vs_regular():
	'''older name:more than 1
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	#general_sc_dict = dict() #key:thr; val: dict of dif
	general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		#sc_dict = {i : 0 for i in range(2,11)}
		fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			#len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#update_d(sc_dict, key, (len_set_cover - len_greedy)/len_greedy)

			#update_d(fsc_dict, key, (len_f_sc - len_set_cover)/len_set_cover)

			update_d(fsc_dict, key, int(bool((len_f_sc - len_set_cover)))) #1 if there was an update, 0 if not

		#was not used to make the original file:

		#
		fam_sizes = [i for i in range(2,10,2)]
		new_fsc_dict = dict()
		for size in fam_sizes:
			if size != 8:
				new_fsc_dict[(size, size+1)] = (fsc_dict[size] + fsc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_fsc_dict[(size, size+1, size+2)] = (fsc_dict[size] + fsc_dict[size+1] + fsc_dict[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])

		#new_fsc_dict[(10)] =
		#
		#fam_sizes = [i for i in range(2,11)]
		#for size in fam_sizes:
			#sc_dict[size] = sc_dict[size]/num_of_fam[size]
		#	fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		#general_sc_dict[thr] = sc_dict
		general_fsc_dict[thr] = new_fsc_dict
	problematics_f.close()
	pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_regular_per_thr_22_02_coupled_more_than_1', 'wb'))
	#pickle.dump(general_sc_dict, open(path + '/sc_dict_diff_from_greedy_per_thr_07_02', 'wb'))

def fuzzy_vs_regular_more_than_1():
	'''older name:more than 1
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	#general_sc_dict = dict() #key:thr; val: dict of dif
	general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		#sc_dict = {i : 0 for i in range(2,11)}
		fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			#len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#update_d(sc_dict, key, (len_set_cover - len_greedy)/len_greedy)

			#update_d(fsc_dict, key, (len_f_sc - len_set_cover)/len_set_cover)
			if (len_f_sc - len_set_cover) >1:
				print('fuzzy better than int by more than 1')
				exit()
				update_d(fsc_dict, key, 1) #1 if there was an update, 0 if not

		#was not used to make the original file:

		#
		fam_sizes = [i for i in range(2,10,2)]
		new_fsc_dict = dict()
		for size in fam_sizes:
			if size != 8:
				new_fsc_dict[(size, size+1)] = (fsc_dict[size] + fsc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_fsc_dict[(size, size+1, size+2)] = (fsc_dict[size] + fsc_dict[size+1] + fsc_dict[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])

		#new_fsc_dict[(10)] =
		#
		#fam_sizes = [i for i in range(2,11)]
		#for size in fam_sizes:
			#sc_dict[size] = sc_dict[size]/num_of_fam[size]
		#	fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		#general_sc_dict[thr] = sc_dict
		general_fsc_dict[thr] = new_fsc_dict
	problematics_f.close()
	pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_regular_per_thr_22_02_coupled_more_than_1_for_real_more_then_1', 'wb'))
	#pickle.dump(general_sc_dict, open(path + '/sc_dict_diff_from_greedy_per_thr_07_02', 'wb'))


def statistics_versatile_thr_greedy_restricted_sc():
	'''CFD for now
	for each family: it's size, the score by restricted set cover and greedy set cover - for 2,3,4 and 5 sgRNAs. it will be the basis for forther analysis.'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	#general_sc_dict = dict() #key:thr; val: dict of dif
	#general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	num_of_fam = count_fam_sizes() #dictionary with the sizes of the family
	full_info_lst = list()
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		#sc_dict = {i : 0 for i in range(2,11)}
		#fsc_dict = {i : 0 for i in range(2,11)}
	#	num_of_fam = {i : 0 for i in range(2,11)}
		for dir in os.listdir(datapath):
			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			#update_d(num_of_fam, num_of_genes, 1)
			##

			greedy = set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb')))
			greedy2, greedy3, greedy4, greedy5 = calculate_score(greedy[:min(2,len(greedy))], genes_lst), calculate_score(greedy[:min(3,len(greedy))], genes_lst), calculate_score(greedy[:min(4,len(greedy))], genes_lst), calculate_score(greedy[:min(5,len(greedy))], genes_lst)
			#calculate here the score for each sub-group of size 2,3,4 or 5. the result from the call is a tuple: (best_score, group). score is avg number of genes that won't get cover.
			cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			if f_sc:
				update_d(fst_dict, key, len(f_sc))
			if set_cover:
				update_d(sc_dict, key, len(set_cover))
			update_d(resc_dict2, key, cover2[0]) #covers_i[0] is 'best score'
			update_d(resc_dict3, key, cover3[0])
			update_d(resc_dict4, key, cover4[0])
			update_d(resc_dict5, key, cover5[0])
			full_info_lst.append([dir, num_of_genes, greedy2, greedy3, greedy4, greedy5, cover2, cover3, cover4, cover5])

			#full_info_dict[dir] = [num_of_genes, greedy2, greedy3, greedy4, greedy5, resc_dict2, resc_dict3, resc_dict4, resc_dict5]

			###

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#
			update_d(sc_dict, key, (len_set_cover - len_greedy)/len_greedy)
			update_d(fsc_dict, key, (len_f_sc - len_greedy)/len_greedy)
		full_info_lst = pd.DataFrame(full_info_lst, columns=['fam_name', 'num_of_genes', 'greedy2', 'greedy3', 'greedy4', 'greedy5', 'cover2', 'cover3', 'cover4', 'cover5'])
		full_info_lst.to_pickle(path + '\full_info_pandas_07_02.p')
		return
		#pickle.save(full_info_lst)
		#was not used to make the original file:
		fam_sizes = [i for i in range(2,11)]
		for size in fam_sizes:
			sc_dict[size] = sc_dict[size]/num_of_fam[size]
			fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		general_sc_dict[thr] = sc_dict
		general_fsc_dict[thr] = fsc_dict
	problematics_f.close()
	pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_greedy_per_thr_06_02', 'wb'))
	pickle.dump(general_sc_dict, open(path + '/sc_dict_diff_from_greedy_per_thr_06_02', 'wb'))


def statistics_versatile_thr_greedy_restricted_sc_fraction():
	'''not finished yet - 0403 CFD for now
	for each family: it's size, the score by restricted set cover and greedy set cover - for 2,3,4 and 5 sgRNAs. it will be the basis for forther analysis.'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	num_of_fam = count_fam_sizes() #dictionary with the sizes of the family
	general_d2, general_d3, general_d4, general_d5 = dict(), dict(), dict(), dict()
	for thr in thrs:
		d2, d3, d4, d5 = {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}
		for dir in os.listdir(datapath):
			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			genes_lst = pickle.load((open(genes_path, "rb")))
			#candidates_lst = pickle.load(open(candidates_path,'rb')
			num_of_genes = len(genes_lst)
			greedy = set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb')))
			#fuzzy_greedy = cover2.gready_cover(,genes_lst )
			greedy2, greedy3, greedy4, greedy5 = calculate_score(greedy[:min(2,len(greedy))], genes_lst), calculate_score(greedy[:min(3,len(greedy))], genes_lst), calculate_score(greedy[:min(4,len(greedy))], genes_lst), calculate_score(greedy[:min(5,len(greedy))], genes_lst)
			#calculate here the score for each sub-group of size 2,3,4 or 5. the result from the call is a tuple: (best_score, group). score is avg number of genes that won't get cover.
			cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.

			update_d(d2, key, int(bool(cover2[0] - greedy2 > 0))) #covers_i[0] is 'best score'
			update_d(d3, key, int(bool(cover3[0] - greedy3 > 0)))
			update_d(d4, key, int(bool(cover4[0] - greedy4 > 0)))
			update_d(d5, key, int(bool(cover5[0] - greedy5 > 0)))


		fam_sizes = [i for i in range(2,10,2)]
		new_d2, new_d3, new_d4, new_d5 = dict(), dict(), dict(), dict()
		for size in fam_sizes:
			if size != 8:
				new_d2[(size, size+1)] = (d2[size] + d2[size+1])/(num_of_fam[size] + num_of_fam[size+1])
				new_d3[(size, size+1)] = (d3[size] + d3[size+1])/(num_of_fam[size] + num_of_fam[size+1])
				new_d4[(size, size+1)] = (d4[size] + d4[size+1])/(num_of_fam[size] + num_of_fam[size+1])
				new_d5[(size, size+1)] = (d5[size] + d5[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_d2[(size, size+1, size+2)] = (d2[size] + d2[size+1] + d2[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])
				new_d3[(size, size+1, size+2)] = (d3[size] + d3[size+1] + d3[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])
				new_d4[(size, size+1, size+2)] = (d4[size] + d4[size+1] + d4[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])
				new_d5[(size, size+1, size+2)] = (d5[size] + d5[size+1] + d5[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])


		general_d2[thr] = new_d2
		general_d3[thr] = new_d3
		general_d4[thr] = new_d4
		general_d5[thr] = new_d5


	problematics_f.close()
	pickle.dump(general_d2, open(path + '/d2_restricted_fsc_NIP_vs_greedy_per_thr_05_03', 'wb'))
	pickle.dump(general_d3, open(path + '/d3_restricted_fsc_MIP_vs_greedy_per_thr_05_03', 'wb'))
	pickle.dump(general_d4, open(path + '/d4_restricted_fsc_MIP_vs_greedy_per_thr_05_03', 'wb'))
	pickle.dump(general_d5, open(path + '/d5_restricted_fsc_MIP_vs_greedy_per_thr_05_03', 'wb'))


def statistics_versatile_thr_restricted_sc_fraction_greedy_fsc_vs_MIP_fsc():
	'''not finished yet - 0603 CFD for now
	for each family: it's size, the score by restricted set cover and greedy set cover - for 2,3,4 and 5 sgRNAs. it will be the basis for forther analysis.'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	num_of_fam = count_fam_sizes() #dictionary with the sizes of the family
	general_d2, general_d3, general_d4, general_d5 = dict(), dict(), dict(), dict()
	for thr in thrs:
		d2, d3, d4, d5 = {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}
		num_of_fam2, num_of_fam3, num_of_fam4, num_of_fam5 = {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}, {i : 0 for i in range(2,11)}
		for dir in os.listdir(datapath):
			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			genes_lst = pickle.load((open(genes_path, "rb")))
			candidates_lst = pickle.load(open(candidates_path,'rb'))
			num_of_genes = len(genes_lst)
			#greedy = set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb')))
			greedy2, greedy3, greedy4, greedy5 = 2 - Covers2.gready_cover(candidates_lst,genes_lst,3)[0] ,3 - Covers2.gready_cover(candidates_lst,genes_lst,3)[0], 4 - Covers2.gready_cover(candidates_lst,genes_lst,4)[0], 5 - Covers2.gready_cover(candidates_lst,genes_lst, 5)[0]
			#greedy2, greedy3, greedy4, greedy5 = calculate_score(greedy[:min(2,len(greedy))], genes_lst), calculate_score(greedy[:min(3,len(greedy))], genes_lst), calculate_score(greedy[:min(4,len(greedy))], genes_lst), calculate_score(greedy[:min(5,len(greedy))], genes_lst)
			#calculate here the score for each sub-group of size 2,3,4 or 5. the result from the call is a tuple: (best_score, group). score is avg number of genes that won't get cover.
			cover2, cover3, cover4, cover5 = 2 - Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2', genes_sg_path, sg_genes_path)[0], 3 - Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path)[0], 4 - Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path)[0], 5 - Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)[0]
			key = num_of_genes #think that kay sepose to be num of gene family.

			if greedy2 < 1.99:
				update_d(d2, key, int(bool(cover2 - greedy2 > 0))) #covers_i[0] is 'best score'
				update_d(num_of_fam2, key, 1)
			if greedy3 < 2.99:
				update_d(d3, key, int(bool(cover3 - greedy3 > 0))) #covers_i[0] is 'best score'
				update_d(num_of_fam3, key, 1)
			if greedy4 < 3.99:
				update_d(d4, key, int(bool(cover4 - greedy4 > 0))) #covers_i[0] is 'best score'
				update_d(num_of_fam4, key, 1)
			if greedy5 < 4.99:
				update_d(d5, key, int(bool(cover5 - greedy5 > 0))) #covers_i[0] is 'best score'
				update_d(num_of_fam5, key, 1)
#			update_d(d3, key, int(bool(cover3[0] - greedy3 > 0)))
#			update_d(d4, key, int(bool(cover4[0] - greedy4 > 0)))
#			update_d(d5, key, int(bool(cover5[0] - greedy5 > 0)))


		fam_sizes = [i for i in range(2,10,2)]
		new_d2, new_d3, new_d4, new_d5 =  dict(), dict(), dict(), dict()
		for size in fam_sizes:
			if size != 8:
				new_d2[(size, size+1)] = (d2[size] + d2[size+1])/(num_of_fam2[size] + num_of_fam2[size+1])
				new_d3[(size, size+1)] = (d3[size] + d3[size+1])/(num_of_fam3[size] + num_of_fam3[size+1])
				new_d4[(size, size+1)] = (d4[size] + d4[size+1])/(num_of_fam4[size] + num_of_fam4[size+1])
				new_d5[(size, size+1)] = (d5[size] + d5[size+1])/(num_of_fam5[size] + num_of_fam5[size+1])
			else:
				new_d2[(size, size+1, size+2)] = (d2[size] + d2[size+1] + d2[size+2])/(num_of_fam2[size] + num_of_fam2[size+1] + num_of_fam2[size+2])
				new_d3[(size, size+1, size+2)] = (d3[size] + d3[size+1] + d3[size+2])/(num_of_fam3[size] + num_of_fam3[size+1] + num_of_fam3[size+2])
				new_d4[(size, size+1, size+2)] = (d4[size] + d4[size+1] + d4[size+2])/(num_of_fam4[size] + num_of_fam4[size+1] + num_of_fam4[size+2])
				new_d5[(size, size+1, size+2)] = (d5[size] + d5[size+1] + d5[size+2])/(num_of_fam5[size] + num_of_fam5[size+1] + num_of_fam5[size+2])


		general_d2[thr] = new_d2
		general_d3[thr] = new_d3
		general_d4[thr] = new_d4
		general_d5[thr] = new_d5


	problematics_f.close()
	pickle.dump(general_d2, open(path + '/d2_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03', 'wb'))
	pickle.dump(general_d3, open(path + '/d3_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03', 'wb'))
	pickle.dump(general_d4, open(path + '/d4_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03', 'wb'))
	pickle.dump(general_d5, open(path + '/d5_restricted_fsc_MIP_vs_greedy__fsc_per_thr_06_03', 'wb'))


def calculate_score(lst_of_candidates, genes_lst):
	'''
	:param lst_of_candidates: list of candidates of the current group
	:param genes_lst: list of the genes of the family
	:return: the objective for this group: sigma pi 1-phi(sg_j, gene)*Xj, i.e, expected sum of genes that wont be cleaves - lower is better
	'''
	#cdef int score, gene_score
	print('lst_of_candidates', lst_of_candidates)
	score = 0
	#print(len(lst_of_candidates), lst)
	for gene in genes_lst:
		gene_score = 1
		for candidate in lst_of_candidates:
			if gene in candidate.genes_score_dict:
				gene_score *= (1- candidate.genes_score_dict[gene])
		#gene_score = 1- gene_score
		score += gene_score
	return score


def count_fam_sizes(datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"):
	#sizes = [i for i in range(2,11)]
	num_of_fam = {i: 0 for i in range(2,11)}#key:thr; val: dict of dif
	for dir in os.listdir(datapath):
		genes_path = "/".join([datapath, dir, "genesNames.p"])
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		update_d(num_of_fam, num_of_genes, 1)
	return num_of_fam

def devide_by_num_of_fam(num_of_fam_dict, res_dict):
	fam_sizes = [i for i in range(2,11)]
	thrs = np.append(np.arange(10)/10, 0.99)
	for t in thrs:
		for size in fam_sizes:
			res_dict[t][size] = res_dict[t][size]/num_of_fam_dict[size]

def from_vers_thr_d_to_table_fuzzy_vs_regular(name, outpath):
	fsc_d = pickle.load(open(outpath + name, 'rb'))
	header = ['fam size', 'strategy',0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
	data = list()
	#fam_size = [i for i in range(2,11)]
	thrs = np.append(np.arange(10)/10, 0.99)
	fam_size = sorted(list(fsc_d[thrs[1]].keys()), key = lambda x: x[0])
	#print(fam_size)
	for size in fam_size:
		fsc_line = [size, 'fsc']
		for t in thrs:
			fsc_line.append(fsc_d[t][size])
		data.append(fsc_line)
	df = pd.DataFrame(data, columns=header)
	df.to_csv(outpath + name + '.csv') #old name by mistake: more than 1

	#df.to_csv(outpath +'/how_many_int_better_greedy_copuled_19_02.csv')

	#df.to_csv(outpath +'/fsc_vs_regular_fraction_22_02.csv') #old name by mistake: more than 1
	#df.to_csv(outpath +'/fsc_vs_regular_more_than_one_04_03.csv') #old name by mistake: more than 1


def from_vers_thr_d_to_table(sc_d_path,fsc_d_path, outpath):
	sc_d = pickle.load(open(sc_d_path, 'rb'))
	fsc_d = pickle.load(open(fsc_d_path, 'rb'))
	num_of_fam = count_fam_sizes()
	#devide_by_num_of_fam(num_of_fam, sc_d)
	#devide_by_num_of_fam(num_of_fam, fsc_d)
	#print(sc_d)
	header = ['fam size', 'strategy',0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
	data = list()

	fam_size = [i for i in range(2,11)]
	#thrs = np.append(np.arange(0.0, 1.0, 0.1), 0.99)
	#thrs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]
	thrs = np.append(np.arange(10)/10, 0.99)
	#print(thrs)
	for size in fam_size:
		sc_line = [size, 'sc']
		fsc_line = [size, 'fsc']
		for t in thrs:
			#sc_num_of_sg = sc_d[t]
			#sc_num_of_sg = sc_num_of_sg[size]
			#fsc_num_of_sg = fsc_d[t][size]
			sc_line.append(sc_d[t][size])
			fsc_line.append(fsc_d[t][size])
		data.append(sc_line)
		data.append(fsc_line)
	df = pd.DataFrame(data, columns=header)
	df.to_csv(outpath +'/comare_covers_08_02.csv')



	#for thr, thr_d in sc_d.items():
	#	for fam_size in thr_d.keys():
	#		sc_num_of_sg = thr_d[fam_size]




			
			
			
	###################
			#find the statistics per this
	#	dir_path = path + str(thr)
	#	if not os.path.exists(dir_path):
	#		os.makedirs(dir_path)
	#	tables_for_greedy_set_cover(datapath, dir_path, thr)
		#tables_of_covers(datapath, dir_path, thr)

def how_many_int_outpreformed_greedy_more_than_1():
	'''CFD for now
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	general_sc_dict = dict() #key:thr; val: dict of dif
	#general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		sc_dict = {i : 0 for i in range(2,11)}
		#fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			#len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#
			if len_set_cover - len_greedy > 1:
				print(len_set_cover - len_greedy)
				exit()
			if len_set_cover - len_greedy > 0:
				update_d(sc_dict, key, int(bool(len_set_cover - len_greedy - 1))) #1 if there was an update, 0 if not
			#update_d(fsc_dict, key, (len_f_sc - len_greedy)/len_greedy)
		#was not used to make the original file:
		fam_sizes = [i for i in range(2,10,2)]
		new_sc_dict = dict()

		for size in fam_sizes:
			if size != 8:
				new_sc_dict[(size, size+1)] = (sc_dict[size] + sc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_sc_dict[(size, size+1, size+2)] = (sc_dict[size] + sc_dict[size+1] + sc_dict[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])

		#for size in fam_sizes:
		#	new_sc_dict[(size, size+1)] = (sc_dict[size] + sc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
		#	fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		general_sc_dict[thr] = new_sc_dict
		#general_fsc_dict[thr] = fsc_dict
	problematics_f.close()
	#pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_greedy_per_thr_07_02', 'wb'))
	pickle.dump(general_sc_dict, open(path + '/how_many_int_better_greedy_more_then_1_copuled_20_02', 'wb'))

def how_many_int_outpreformed_greedy():
	'''CFD for now
	compute the relative diference in number of sgRNA being used. yields avg and SD for each fam size, for each thr'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_05_02"#"/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"#
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr"
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	general_sc_dict = dict() #key:thr; val: dict of dif
	#general_fsc_dict = dict() #key:thr; val: dict of dif
	problematics_f = open(path+'/problematics_04_02.txt', 'w')
	for thr in thrs:
		#greedy_dict = {i : 0 for i in range(2,11)}
		sc_dict = {i : 0 for i in range(2,11)}
		#fsc_dict = {i : 0 for i in range(2,11)}
		num_of_fam = {i : 0 for i in range(2,11)}

		#print('here0')
		for dir in os.listdir(datapath):
			#print('here1')

			candidates_path = "/".join([datapath, dir, "res_in_lst.p"])
			genes_path = "/".join([datapath, dir, "genesNames.p"])
			genes_sg_path = "/".join([datapath, dir, "genes_sg_dict.p"])
			sg_genes_path = "/".join([datapath, dir, "sg_genes_dict.p"])
			if not os.path.isfile(candidates_path):
				problematics_f.write(dir + ": no cadideates lst\n")
				continue
			num_of_genes = len(pickle.load((open(genes_path, "rb"))))
			#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
			#if not best_thr_candidate:
			#	problematics_f.write(dir + ": empty candidates lst\n")
			#	continue
			#if best_amount == num_of_genes:
				#update_d(full_cover,num_of_genes,1)

			update_d(num_of_fam, num_of_genes, 1)

			len_set_cover = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "SC", genes_sg_path, sg_genes_path))
			len_greedy = len(set_cover_greedy.find_set_cover(pickle.load(open(candidates_path,'rb')), pickle.load(open(sg_genes_path,'rb')), thr, pickle.load(open(genes_sg_path,'rb'))))
			#len_f_sc = len(Covers2.call_CplexCovers(candidates_path, genes_path, thr, "F_SC", genes_sg_path, sg_genes_path))
			#cover2, cover3, cover4, cover5 = Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_2',genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_3', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_4', genes_sg_path, sg_genes_path), Covers2.call_CplexCovers(candidates_path, genes_path, 0.99, 'BC_5', genes_sg_path, sg_genes_path)
			key = num_of_genes #think that kay sepose to be num of gene family.
			#not sure the next 3 are needed
			#update_d_SD(greedy_dict, key, greedy)
			#update_d_SD(sc_dict, key, greedy)
			#update_d_SD(fst_dict, key, greedy)
			#
			update_d(sc_dict, key, int(bool((len_set_cover - len_greedy)))) #1 if there was an update, 0 if not
			#update_d(fsc_dict, key, (len_f_sc - len_greedy)/len_greedy)
		#was not used to make the original file:
		fam_sizes = [i for i in range(2,10,2)]
		new_sc_dict = dict()

		for size in fam_sizes:
			if size != 8:
				new_sc_dict[(size, size+1)] = (sc_dict[size] + sc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
			else:
				new_sc_dict[(size, size+1, size+2)] = (sc_dict[size] + sc_dict[size+1] + sc_dict[size+2])/(num_of_fam[size] + num_of_fam[size+1] + num_of_fam[size+2])

		#for size in fam_sizes:
		#	new_sc_dict[(size, size+1)] = (sc_dict[size] + sc_dict[size+1])/(num_of_fam[size] + num_of_fam[size+1])
		#	fsc_dict[size] = fsc_dict[size]/num_of_fam[size]
		#
		general_sc_dict[thr] = new_sc_dict
		#general_fsc_dict[thr] = fsc_dict
	problematics_f.close()
	#pickle.dump(general_fsc_dict, open(path + '/fsc_dict_diff_from_greedy_per_thr_07_02', 'wb'))
	pickle.dump(general_sc_dict, open(path + '/how_many_int_better_greedy_copuled_19_02', 'wb'))

#def merege_fzm_sizes(input_dict):
#	'''
#	dict_res: dict of dicts.
#	key: thr, val: dict - key:fam size, val - value
#	'''
	#num_of_fam = count_fam_sizes(datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11")
	#for thr in input_dict:
	#	for fam_size in input_dict[thr]:

	#return

def update_d(d, key, val):
	if d[key] == 0:
		d[key] = val
	else:
		d[key] = d[key] + val
		
def update_d_27_1(d, key, val):
	if key in d:
		d[key] = d[key] + val
	else:
		d[key] = val
		
def update_d_SD(d, key, val):
	'''in order to compute the SD, all of the values are saved. At the end, SD is computed
	val is a list o values'''
	if key in d:
		d[key] = d[key] + [val]
	else:
		d[key] = [val]

def how_many_full_cover(SD_above_thr_dict):
	for key, val in SD_above_thr_dict.items():
		full = 0
		for v in val:
			full += (lambda x: x == key)(v)
		full = full/len(val)
		print("fam size: " + str(key) + "full cover: " + str(full))
	
def find_thr_best(path, outpath, thr):
	'''find the best sgRNA by the maximising number of genes being cleaved above thr, and summeraize the data to a table'''
	full_cover = {}
	num_of_fam = {}
	cut_expectation = {}
	SD_dict = dict()
	f2_path = "/".join([outpath, "thr_best_analysis_04_11.txt"])
	problematics_f = open("/".join([outpath, "problematics_thr_best_analysis_04_11.txt"]),'w')
	for dir in os.listdir(path):
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		#best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		#if not best_thr_candidate:
		#	problematics_f.write(dir + ": empty candidates lst\n")
		#	continue
		#if best_amount == num_of_genes:
			#update_d(full_cover,num_of_genes,1)
		update_d(full_cover,num_of_genes,best_amount == num_of_genes)
		update_d(num_of_fam, num_of_genes, 1)
		update_d(cut_expectation, num_of_genes, best_thr_candidate.cut_expectation)
		update_d_SD(SD_dict, num_of_genes, best_thr_candidate.cut_expectation)
			##write to to file
	problematics_f.close()
	with open(f2_path, 'w') as csvfile:
			#spamwriter = csv.writer(csvfile, delimiter=' ',
			#                quotechar='|', quoting=csv.QUOTE_MINIMAL)
			#spamwriter.writerow(["fam size", "number of fam", "fraction of full cover by the first sgRNA"])
			#for i in range(2,11): #the sizes of the families
			#    spamwriter.writerow([str(i), str(num_of_fam[i]), str(full_cover[i]/num_of_fam[i])])
			csvfile.write("fam size;number of fam;fraction of full cover by the first sgRNA;avg_cut_expectation\n")
			for i in range(2,11): #the sizes of the families
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				if i in num_of_fam:
					csvfile.write(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i])+ ";"+ str(cut_expectation[i]/num_of_fam[i]) + "\n")
					#print("SD of size " + str(i) +": " +  str(numpy.std(numpy.array(SD_dict[i]))))

def find_fraction_of_cover(path, outpath, thr):
	'''find how many s1 candidates are predicted to cleave all of the gene family'''
	full_cover = {}
	num_of_fam = {}
	cut_expectation = {}
	f2_path = "/".join([outpath, "num_of_covered_04_11.txt"])
	problematics_f = open("/".join([outpath, "problematics_num_of_covers_04_11.txt"]),'w')
	for dir in os.listdir(path):
		#print(dir)
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		candidate, amount = Covers2.find_fraction(candidates_path, thr)
		if not candidate:
			problematics_f.write(dir + ": empty candidates lst\n")
			continue
#		if amount == num_of_genes:
#			update_d(full_cover,num_of_genes,1)
		update_d(full_cover,num_of_genes, amount == num_of_genes)
		update_d(num_of_fam, num_of_genes, 1)
		update_d(cut_expectation, num_of_genes, candidate.cut_expectation)
			##write to to file
	problematics_f.close()
	with open(f2_path, 'w') as csvfile:
			#spamwriter = csv.writer(csvfile, delimiter=' ',
			#                quotechar='|', quoting=csv.QUOTE_MINIMAL)
			#spamwriter.writerow(["fam size", "number of fam", "fraction of full cover by the first sgRNA"])
			#for i in range(2,11): #the sizes of the families
			#    spamwriter.writerow([str(i), str(num_of_fam[i]), str(full_cover[i]/num_of_fam[i])])
			csvfile.write("fam size;number of fam;fraction of full cover by the first sgRNA;avg_cut_expectation\n")
			for i in range(2,11): #the sizes of the families
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				if i in num_of_fam:
					csvfile.write(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i])+ ";"+ str(cut_expectation[i]/num_of_fam[i]) + "\n")

def find_number_of_genes_covered_by_first(path, outpath, thr, s):
	'''find the number of genes cleaved w.p>thr per family, for s1 or s2 candidates'''
	cover_dict = {}
	num_of_fam = {}
	cut_expectation = {}
	SD_above_thr_dict = dict()
	SD_expectation_dict = dict()
	f2_path = "/".join([outpath, "number_of_genes_covered_by_first_04_11.txt"])
	problematics_f = open("/".join([outpath, "problematics_num_of_covers_04_11.txt"]),'w')
	counter = 0
	for dir in os.listdir(path):
		counter += 1
		#print(dir)
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		if s == 's1':
			candidate, amount = Covers2.find_fraction(candidates_path, thr) #how many genes above omega

		elif s == 's2':
			candidate, amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		if not candidate:
			problematics_f.write(dir + ": empty candidates lst\n")
			continue
		#print(dir)
		#print(amount)
		#print(num_of_genes)
		#print('genes list:', pickle.load((open(genes_path, "rb"))))
		#if amount == num_of_genes:
		update_d(cover_dict, num_of_genes, amount)
		update_d(num_of_fam, num_of_genes, 1)
		update_d(cut_expectation, num_of_genes, candidate.cut_expectation)
		update_d_SD(SD_expectation_dict, num_of_genes, candidate.cut_expectation)
		update_d_SD(SD_above_thr_dict, num_of_genes, amount)
		##write to to file
	problematics_f.close()
	print("counter: ", counter)
	with open(f2_path, 'w') as csvfile:
			#spamwriter = csv.writer(csvfile, delimiter=' ',
			#                quotechar='|', quoting=csv.QUOTE_MINIMAL)
			#spamwriter.writerow(["fam size", "number of fam", "fraction of full cover by the first sgRNA"])
			#for i in range(2,11): #the sizes of the families
			#    spamwriter.writerow([str(i), str(num_of_fam[i]), str(full_cover[i]/num_of_fam[i])])
			csvfile.write("fam size;number of fam;fraction of full cover by the first sgRNA;avg_cut_expectation\n")
			#print("fam size;number of fam;fraction of full cover by the first sgRNA;avg_cut_expectation\n")
			how_many_full_cover(SD_above_thr_dict)

			for i in range(2,11): #the sizes of the families
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				print(num_of_fam)#" +";"+ str(cover_dict[i]/num_of_fam[i]) + "(" +str(numpy.std(numpy.array(SD_above_thr_dict[i]))) +");"+ str(cut_expectation[i]/num_of_fam[i]) + "(" + str(numpy.std(numpy.array(SD_expectation_dict[i]))) + ")\n")
				csvfile.write(str(i) +";"+ str(num_of_fam[i]) +";"+ str(cover_dict[i]/num_of_fam[i]) + "(" +str(numpy.std(numpy.array(SD_above_thr_dict[i]))) +");"+ str(cut_expectation[i]/num_of_fam[i]) + "(" + str(numpy.std(numpy.array(SD_expectation_dict[i]))) + ")\n")
				#print("SD expectation of size " + str(i) +": " +  str(numpy.std(numpy.array(SD_expectation_dict[i]))))
				#print("SD above thr of size " + str(i) +": " +  str(numpy.std(numpy.array(SD_above_thr_dict[i]))))
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(cover_dict[i]/num_of_fam[i])+ ";"+ str(cut_expectation[i]/num_of_fam[i]) + "\n")

def data_for_r_barplot(path, outpath, thr):
	covers_array = numpy.zeros((11,11), dtype = int) # = [[0 for i in range(10)] for i in range(10)]
	num_of_fam = numpy.zeros(11)# = [0 for i in range(10)]
	shiran_style_matrix = numpy.zeros((11,11), dtype = numpy.int)
	#cut_expectation = {}
	f2_path = "/".join([outpath, "best_sg_num_of_cleaved_04_11.csv"])
	problematics_f = open("/".join([outpath, "problematics_num_of_covers_426036.txt"]),'w')
	for dir in os.listdir(path):
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
	#	if s == 's1':
	#		candidate, amount = Covers2.find_fraction(candidates_path, thr)
	#	elif s == 's2':
		candidate, amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)

		if not candidate:
		#	problematics_f.write(dir + ": empty candidates lst\n")
			continue
		#if amount == num_of_genes:
		covers_array[amount][num_of_genes] += 1
		num_of_fam[num_of_genes] += 1
	#normalized
	print(covers_array)
	shiran_style_matrix = covers_array[1:9, 2:11]
	#write to a file
	f2 = open(f2_path, 'w')
	f2.write("# cleaved by best sg | family size\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n")
	f2.write(str(covers_array))

	numpy.savetxt("np_savetxt.txt", covers_array, delimiter = '\t',  header="guides per family|family size,2,3,4,5,6,7,8,9,10")
	rows_titles = numpy.transpose(numpy.asmatrix(numpy.array(list(range(1,9)))))
	numpy.savetxt("shiran_style_matrix_09_04.csv", numpy.hstack((rows_titles, shiran_style_matrix)), fmt='%i', delimiter=",", header="guides per family|family size,2,3,4,5,6,7,8,9,10")

		#update_d(cover_dict, num_of_genes, amount)
		#update_d(num_of_fam, num_of_genes, 1)
		#update_d(cut_expectation, num_of_genes, candidate.cut_expectation)
			##write to to file
	#problematics_f.close()
	#with open(f2_path, 'w') as csvfile:
			#spamwriter = csv.writer(csvfile, delimiter=' ',
			#                quotechar='|', quoting=csv.QUOTE_MINIMAL)
			#spamwriter.writerow(["fam size", "number of fam", "fraction of full cover by the first sgRNA"])
			#for i in range(2,11): #the sizes of the families
			#    spamwriter.writerow([str(i), str(num_of_fam[i]), str(full_cover[i]/num_of_fam[i])])
	
	#	csvfile.write("fam size;number of fam;fraction of full cover by the first sgRNA;avg_cut_expectation\n")
	#		for i in range(2,11): #the sizes of the families
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
				#print(str(i) +";"+ str(num_of_fam[i])+";"+ str(full_cover[i]/num_of_fam[i]))
	#			if i in num_of_fam:
	#				csvfile.write(str(i) +";"+ str(num_of_fam[i])+";"+ str(cover_dict[i]/num_of_fam[i])+ ";"+ str(cut_expectation[i]/num_of_fam[i]) + "\n")
	f2.close()
	

def data_for_r_barplot_only_amount(path, outpath, thr):
	covers_array = numpy.zeros((11,11), dtype = int) # = [[0 for i in range(10)] for i in range(10)]
	num_of_fam = numpy.zeros(11)# = [0 for i in range(10)]
	shiran_style_matrix = numpy.zeros((11,11), dtype = numpy.int)
	#cut_expectation = {}
	f2_path = "/".join([outpath, "best_sg_num_of_cleaved_04_11.csv"])
	#problematics_f = open("/".join([outpath, "problematics_num_of_covers_426036.txt"]),'w')
	for dir in os.listdir(path):
		#print(dir)
		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		genes_path = "/".join([path, dir, "genesNames.p"])
		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		num_of_genes = len(pickle.load((open(genes_path, "rb"))))
	#	if s == 's1':
	#		candidate, amount = Covers2.find_fraction(candidates_path, thr)
	#	elif s == 's2':
		#candidate, amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		amount = Covers2.test_thr_best_only_amount(candidates_path)
		#if not candidate:
		if not amount:
		#	problematics_f.write(dir + ": empty candidates lst\n")
			continue
		#if amount == num_of_genes:
		covers_array[amount][num_of_genes] += 1
		num_of_fam[num_of_genes] += 1
	#normalized
	shiran_style_matrix = covers_array[1:9, 2:11]
	#write to a file
	f2 = open(f2_path, 'w')
	f2.write("# cleaved by best sg | family size\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n")
	f2.write(str(covers_array))

	numpy.savetxt("np_savetxt.txt", covers_array, delimiter = '\t',  header="guides per family|family size,2,3,4,5,6,7,8,9,10")
	rows_titles = numpy.transpose(numpy.asmatrix(numpy.array(list(range(1,9)))))
	numpy.savetxt("shiran_style_matrix_09_04.csv", numpy.hstack((rows_titles, shiran_style_matrix)), fmt='%i', delimiter=",", header="guides per family|family size,2,3,4,5,6,7,8,9,10")

def how_many_more_then_one(path):
	res, num_of_fam = 0, 0
	for dir in os.listdir(path):
		num_of_fam += 1
		#print(dir)
		candidates_path = "/".join([path, dir, "res_in_lst.p"])

		if not os.path.isfile(candidates_path):
			continue
			
		candidates_lst = pickle.load(open(candidates_path, 'rb'))
		genes_path = "/".join([path, dir, "genesNames.p"])

		for c in candidates_lst:
			for v in c.targets_dict.values():
				if len(v) > 1:
					print(c)
					print('\n\n')
					res += 1
	return res
	
def find_difrances_in_runs(run1path, raun2path):
	return
	
def versatile_thr_covers():
	'''CFD for now'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr" 
	thrs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	#thrs = [0.0, 0.5]
	for thr in thrs:
		print('thr: ', thr)
		dir_path = path + str(thr)
		if not os.path.exists(dir_path):
			os.makedirs(dir_path)
		tables_for_greedy_set_cover(datapath, dir_path, thr)
		#tables_of_covers(datapath, dir_path, thr)

	
def analyse_versatile_thr():
	'''CFD for now'''
	datapath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"
	#path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/greedy_set_cover/versatile_thr" 
	thrs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]#change with the corespondig precentiles
	#thrs = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0]#change with the corespondig precentiles
	thrs = [0.0, 0.5]
	for thr in thrs:
		print('thr: ', thr)
		dir_path = path + str(thr)
		if not os.path.exists(dir_path):
			os.makedirs(dir_path)
		outpath_thr_best = dir_path + "/ThrBestAnalysis/"
		if not os.path.exists(outpath_thr_best):
			os.makedirs(outpath_thr_best)
		outpath_fraction_of_cover = dir_path + "/cutExpectationAnalysis/"
		if not os.path.exists(outpath_fraction_of_cover):
			os.makedirs(outpath_fraction_of_cover)
			#outpath_tables_of_covers = dir_path + "/tables_of_cover/"
			#if not os.path.exists(outpath_tables_of_covers):
			#	os.makedirs(outpath_tables_of_covers)
	
		#tables_of_covers(datapath, outpath_tables_of_covers, thr)
		find_thr_best(datapath, outpath_thr_best, thr)
		#print("first function")
		find_fraction_of_cover(datapath, outpath_fraction_of_cover, thr)
		print("fraction of cover")
		find_number_of_genes_covered_by_first(datapath, outpath_fraction_of_cover, thr, "s1")
		print("thr best")
		find_number_of_genes_covered_by_first(datapath, outpath_thr_best, thr, "s2")
		data_for_r_barplot(datapath, outpath_thr_best, thr)

if __name__ == '__main__':
	#statistics_versatile_thr_restricted_sc_fraction_greedy_fsc_vs_MIP_fsc()
	#exit()
	#statistics_versatile_thr_greedy_restricted_sc_fraction()
	from_vers_thr_d_to_table_fuzzy_vs_regular('d2_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')
	from_vers_thr_d_to_table_fuzzy_vs_regular('d3_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')
	from_vers_thr_d_to_table_fuzzy_vs_regular('d4_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')
	from_vers_thr_d_to_table_fuzzy_vs_regular('d5_restricted_fsc_MIP_vs_greedy_fsc_per_thr_06_03','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')

	exit()
	#test_compare_greedy_regular()
	#exit()
	#diff_versatile_thr_greedy_regular_fuzzy_covers()
	#from_vers_thr_d_to_table('/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/sc_dict_diff_from_greedy_per_thr_07_02', '/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/fsc_dict_diff_from_greedy_per_thr_07_02' ,'/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')
	#how_many_int_outpreformed_greedy()
	#fraction_fuzzy_vs_regular()
	#fuzzy_vs_regular_more_than_1()
	#exit()
	#how_many_int_outpreformed_greedy_more_than_1()
	#exit()
	#from_vers_thr_d_to_table_fuzzy_vs_regular('/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/fsc_dict_diff_from_regular_per_thr_20_02_coupled','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')
	from_vers_thr_d_to_table_fuzzy_vs_regular('/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/fsc_dict_diff_from_regular_per_thr_22_02_coupled_more_than_1','/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/set_cover/versatile_thr/')

	exit()
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_04_11"
	outpath = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/greedy_set_cover"
	thr = 0.45#0.426036	
	tables_for_greedy_set_cover(path, outpath, thr)
	exit()

	#path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6"
	#outpath_thr_best = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_thr_best"
	#outpath_fraction_of_cover = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_fraction_of_cover"
	#path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_BU_V1_6_0609"
	#path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409"
	
	#outpath_thr_best = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/CFD/ThrBestAnalysis_0_43/"
	#outpath_fraction_of_cover = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/CFD/cutExpectationAnalysis_0_43/"
	outpath_tables_of_covers = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/MetaAnalysis/CFD/tables_of_cover_0_45/"

	
	#print(how_many_more_then_one("/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/Sample/SampleData"))
	#thr = 0.33
	thr = 0.45#0.426036
	#thr = 0.66157
	#thr = 0.0
	tables_of_covers(path, outpath_tables_of_covers, thr)
	exit()
	find_thr_best(path, outpath_thr_best, thr)
	#print("first function")
	find_fraction_of_cover(path, outpath_fraction_of_cover, thr)
	print("fraction of cover")
	find_number_of_genes_covered_by_first(path, outpath_fraction_of_cover, thr, "s1")
	print("thr best")
	find_number_of_genes_covered_by_first(path, outpath_thr_best, thr, "s2")
	data_for_r_barplot(path, outpath_thr_best, thr)
	

