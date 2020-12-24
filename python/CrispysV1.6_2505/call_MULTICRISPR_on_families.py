__author__ = 'ItayM5'
import re
import os
import bottemsUpAlgorithm
import CasSites

def make_organism_families_file(org, genefamily_data_path, outpath):
	'''org: SL for solanum l. '''
	f1 = open(genefamily_data_path)
	next(f1)
	fout = open(outpath, "w")
	for line in f1:
		line_as_array = line.split(",")
		line_as_array = line_as_array[0].split(";") + line_as_array[1:]
		#print(line_as_array)
		len_of_array = len(line_as_array)
		i = 1
		while i < len_of_array:
		#for i in range(1, len_of_array):
			if org not in line_as_array[i]:
				line_as_array.pop(i)
				len_of_array -= 1
			else:
				i+=1
		fout.write(str(line_as_array)[1:-1]+ "\n")
	f1.close()
	fout.close()

def find_sgRNAs_for_family(family_name, families_list_path, sgRNA_folder, outpath, members_list = None):
	'''families list path: the file made by "make_organism_families_file", from the fenefamily_data_file'''
	##first: find family members:
	extra = "1.1RNAfile.txt-holeGenomeWithExtras-s2n100pC"
	extraVers2 = "2.1RNAfile.txt-holeGenomeWithExtras-s2n100pC"
	if(family_name):
		f1 = open(families_list_path)
		for line in f1:
			line_as_array = line.split(",")
			if family_name in line_as_array[0]:  ##found our family
				members_list = line_as_array[1:]
				break
	if not members_list:
		print("coldn't find the family in the file\n")
		exit(-1)
	outfile = open(outpath, "w")  ##might won't be used
	res = []  ##list of tuples. each tuple: (geneName, [(sgName, sgSeq)...])
	for gene in members_list:
		sgList = []
		geneName = convert_format(gene)
		sites_path = sgRNA_folder + "\\" + geneName +"\\"+ extra + "\\_sites.txt"
		sites_file = open(sites_path, "r")
		for line in sites_file:
			if line[0] == '>':
				sgName = line  ##need to check if there is a need to remove \n
			else:
				seq = line
				sgList += [(sgName, seq)]  #can remove the inner (), and then, recognise by: unevens are names, evens are sequences
		outfile.write(geneName + "\t"+ str(sgList[1:-1])+ "\n")
		res.append((geneName, [sgList]))
	outfile.close()
	return  res

def convert_format(gene):
	gene = re.sub('G', 'g', gene)
	return re.sub('SL', 'Solyc', gene)

def call_MULTYCRISPER(genes_sg_list):       ##the list: list of tuples. each tuple: (geneName, [(sgName, sgSeq), ...])
	sg_list = []
	names_list = []
	genes_sg_dict = {}  #key: gene name. value: list of sgRNA seq
	sg_genes_dict = {}  #key: sgRNA. values: list of genes name. usually this list will be of size 1. need it for the bottems up algorithm
	for tuple in genes_sg_list: #each tuple represent a gene
		dict_val = []
		for sg in tuple[1]:  #tuple[1] is the list of tuples: [(sgName, sgSeq)...])
			names_list.append(sg[1])  #for now, the name will be the seq.
			sg_list.append((sg[1]))
			dict_val.append(sg[1])
			if sg[1] not in sg_genes_dict:
				sg_genes_dict[sg[1]] = [tuple[0]]
			else:
				sg_genes_dict[sg[1]] = sg_genes_dict[sg[1]].append(tuple[0])

		genes_sg_dict[tuple[0]] = dict_val  ##tuple[0] is the gene name

		##adding sg_genes_dict: key: sgRNA. values: gene name. need it for the bottems up algorithm


	##here, make the call. it might cange to 2 seperated functions


def find_sgRNAs_for_family_V2(family_name, families_list_path, sgRNA_folder, members_list = None):
	'''without the need to call_MULTYCRISPER: wifout folding the lists to dictionaries
	families list path: the file made by "make_organism_families_file", from the fenefamily_data_file'''
	##first: find family members:
	extra = "1.1RNAfile.txt-holeGenomeWithExtras-s2n100pC"
	extraVers2 = "2.1RNAfile.txt-holeGenomeWithExtras-s2n100pC"
	if(family_name):
		f1 = open(families_list_path)
		for line in f1:
			line_as_array = line.split(",")
			if family_name in line_as_array[0]:  ##found our family
				members_list = line_as_array[1:]
				break
	if not members_list:
		print("coldn't find the family in the file\n")
		exit(-1)
	#res = []  ##list of tuples. each tuple: (geneName, [(sgName, sgSeq)...])
	sgList = []
	sgNames = []
	genes_sg_dict = {}
	for gene in members_list:
		sgList = []
		sgNames = []
		dict_val = []
		geneName = convert_format(gene)
		sites_path = sgRNA_folder + "\\" + geneName +"\\"+ extra + "\\_sites.txt"
		sites_file = open(sites_path, "r")
		for line in sites_file:
			if line[0] == '>':
				sgName = line  ##need to check if there is a need to remove \n
			else:
				seq = line
				sgList += [seq]  #can remove the inner (), and then, recognise by: unevens are names, evens are sequences
				sgNames += [seq]  ##can be sgName insted. for now it looks like seq is more effitioant. can be change to just make a copy so this at the end
				dict_val += [seq]
		genes_sg_dict[geneName] = dict_val
	return sgList, sgName, genes_sg_dict

def call_on_real_families(paths, splitword):
	#if not path:
	#path = "D:\Gal\MultiCrisper\Eilon familis\gray one\Solyc05g009500.2.1.txt-format.holeGenomeWithExtars-s2pN\_sites.txt"
	genes_sg_dict = {}
	sg_genes_dict = {}
	sgNames = []
	sgList = []
	for p in paths:
		f = open(p)
		res = []
		for line in f:
			if line[0] != ">":
				res += [line[:-4]]
		f.close()
		gene_name = p.split(".txt")[0]
		gene_name = gene_name.split(splitword)[1]
		genes_sg_dict[gene_name] = res
		for sg in res:
			if sg in sg_genes_dict:
				sg_genes_dict[sg] = sg_genes_dict[sg] + [gene_name]
			else:
				sg_genes_dict[sg] = [gene_name]
			if sg not in sgNames:
				sgNames.append(sg)
				sgList.append(sg)
	print(bottemsUpAlgorithm.call_it_all(sgList, sgNames, sg_genes_dict))

def call_E_gray_family():
	dirp = "D:\Gal\MultiCrisper\Eilon familis\gray one"
	paths = [dirp + "\Solyc05g009500.2.1.txt-format.holeGenomeWithExtars-s2pN\_sites.txt", dirp + "\Solyc06g005070.1.1.txt-format.holeGenomeWithExtars-s2pN\_sites.txt"]
	call_on_real_families(paths)

def call_Meirav_family():
	dirp = "D:\Gal\MultiCrisper\Adi families\meirav"
	paths = []
	for dir in os.listdir(dirp):
		if "format" in dir:
			 paths.append(dirp + "\\" +dir+ "\\_sites.txt" )
	call_on_real_families(paths, "meirav\\")

def call_E_purple_family():
	dirp = "D:\Gal\MultiCrisper\Eilon familis\dark_purple"
	paths = []
	for dir in os.listdir(dirp):
		if "format" in dir:
			 paths.append(dirp + "\\" +dir+ "\\_sites.txt" )
		call_on_real_families(paths, "dark_purple\\")

def call_on_given_family(dirp, splitword):
	paths = []
	for dir in os.listdir(dirp):
		if "format" in dir:
			 paths.append(dirp + "\\" +dir+ "\\_sites.txt" )
	call_on_real_families(paths, splitword)

def call_on_E_redundant_purples():
	dirp = "D:\Gal\MultiCrisper\Eilon familis\dark_purple"
	paths = [dirp + "\Solyc01g091180.2.1RNAfile.txt-format.holeGenomeWithExtars-s2pN\_sites.txt", dirp + "\Solyc09g011400.1.1RNAfile.txt-format.holeGenomeWithExtars-s2pN\_sites.txt"]
	call_on_real_families(paths,"dark_purple\\")



	##now, printing the part for the dict:
def call_using_CasSites(dirp, Omega = 0.11, min_length = 20, max_length = 20,start_with_G = False, on_redundant = False, redundant_genes = []):
	'''in dirp will be a list files. In each there will be a sequences in FASTA format'''
	#genes_list = []
	genes_sg_dict = {}
	sg_genes_dict = {}
	sgNames = []
	sgList = []
	for p in os.listdir(dirp):
		#if p[-11:-1]== "RNAfile.tx":  ##g gene file
		if "RNAfile.tx" in p:
			gene_name = p.split(".txt")[0]
			#print(gene_name)
			if(on_redundant):
				if gene_name not in redundant_genes:
					continue
			f = open(dirp + "\\"+ p)
			next(f)
			gene = f.read()
			gene.replace('/n', '')
			#oledr version:
			#for line in f: #only 1 line left
				#genes_list.append(gene_name,line)
			#	gene = line.rstrip()
			f.close()
			genes_sg_dict[gene_name] = CasSites.get_sites(gene, min_length, max_length,start_with_G)
			for sg in genes_sg_dict[gene_name]:
				if sg in sg_genes_dict:
					sg_genes_dict[sg] = sg_genes_dict[sg] + [gene_name]
				else:
					sg_genes_dict[sg] = [gene_name]
				if sg not in sgNames:
					sgNames.append(sg)
					sgList.append(sg)

	return(bottemsUpAlgorithm.call_it_all(sgList, sgNames, sg_genes_dict, Omega))




if __name__ == "__main__":
	'''
	dirp_in = "D:\Gal\MultiCrisper"
	dirp_out = dirp_in + "\output_files"
	org = 'SL'
	genefamily_data_path = dirp_in + "\genefamily_data.hom.csv"
	outpath = dirp_out + "\SLY_genefamily_data.hom.txt"

	make_organism_families_file(org, genefamily_data_path, outpath)
	'''

	#call_E_gray_family()
#	call_Meirav_family()
	#call_on_E_redundant_purples()
	#dirp = "D:\Gal\MCdata\Eilon familis\\arabidopsis\\red\onlyRNA"
	dirp = "D:\Gal\MCdata\Eilon familis\\arabidopsis\\yellow\onlyRNA"
	print(call_using_CasSites(dirp))
	#call_using_CasSites(dirp, True, ['Solyc12g089230.1.1RNAfile', 'Solyc08g061010.2.1RNAfile'])
