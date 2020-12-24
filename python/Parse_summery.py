import analyse_covers
import os
import Covers2
import pickle
import Candidate


def parse_summery(inpath):
	fam_bestSg_dict = dict()
	with open(inpath,'r') as f:
		next(f)
		for line in f:
			l_a_a = line.split(',')
			if line[0][0] != 'H':
				continue
			fam_bestSg_dict[l_a_a[0]] = int(l_a_a[3])
	return fam_bestSg_dict


def best_candidate_from_files(path, thr = 0.66):
	res = {}
	#num_of_fam = {}
	#cut_expectation = {}
	#SD_dict = dict()
	#f2_path = "/".join([outpath, "thr_best_analysis_0_66157.txt"])
	#problematics_f = open("/".join([outpath, "problematics_thr_best_analysis_0_66157.txt"]),'w')
	for dir in os.listdir(path):
		#print(dir)

		candidates_path = "/".join([path, dir, "res_in_lst.p"])
		#genes_path = "/".join([path, dir, "genesNames.p"])

		if not os.path.isfile(candidates_path):
			problematics_f.write(dir + ": no cadideates lst\n")
			continue
		#num_of_genes = len(pickle.load((open(genes_path, "rb"))))
		best_thr_candidate, best_amount, best_cut_them_all = Covers2.test_thr_best(candidates_path, thr)
		if not best_thr_candidate:
			continue
		#if best_amount == num_of_genes:
			#update_d(full_cover,num_of_genes,1)
		res[dir] = best_amount
	return res

def similarity_in_dict(res1,res2, output):
	f = open(output, 'w')
	f.write("old version;new version\n")
	for famname in res1.keys():
		if not famname in res2:
			print(famname)
		elif res1[famname] != res2[famname]:
			f.write(famname + ";" + str(res1[famname]) +";"+ str(res2[famname])+"\n")

			
def show_candidate(path, candidate_name):
	candidates_path = "/".join([path, dir, "res_in_lst.p"])
	c_lst = pickle.load((open(candidates_path, "rb")))
	for c in c_lst:
		if c.name == candidate_name:
			print(cadideates)
			return
	
	

if __name__ == "__main__":
	old_run_file = "summery_CDF_score_Omega_0.66157_new_stopping_cond_11_genes08112016.txt"
	new_run_path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_6_0409"
	outpath = "old_new_diff.txt"
	d1, d2 = parse_summery(old_run_file), best_candidate_from_files(new_run_path)
	pickle.dump(d1, open("old_rin_file.p","wb"))
	pickle.dump(d2, open("new_rin_file.p","wb"))
	similarity_in_dict(d1, d2, outpath)
