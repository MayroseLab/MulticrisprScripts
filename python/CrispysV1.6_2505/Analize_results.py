__author__ = 'ItayM5'
import os

def parse_file(in_file):
	'''fam_size, #guides, #genes per best guide, AVG(genes per guide) #mm of best guide, AVG[#mm per guide-target]'''
	f1 = open(in_file, 'r')
	first_SG = -1
	lines = f1.readlines()
	block_size = 9
	#first guide info
	num_of_gene_cleaved_by_first = len(lines[5].split("  ")) # asldfkjalsdkfjl;askffkl
	#avg num of mm
	AVG_of_mm = 0
	fam_memmbers = set()
	num_of_guides = 0
	i = 0
	while (i+1)*block_size <= len(lines):
	# num mm of first
	#line_7 = lines[7]
		if i == 0:
			AVG_num_of_mm_first_guide = (lines[i*block_size + 7].count(':'))/num_of_gene_cleaved_by_first
		cleaved_genes = lines[i*block_size + 5].split("  ")
		num_of_gene_cleaved = len(cleaved_genes)
		AVG_of_mm += (lines[i*block_size + 7].count(':'))/num_of_gene_cleaved_by_first
		for gene in cleaved_genes:
			if "SL" in gene: #sanity cheack
				fam_memmbers.add(gene.strip())

		i += 1

	#count num of sg
	famsize = len(set)
	num_of_guides = i
	f1.close()
	return([famsize, num_of_guides, num_of_gene_cleaved_by_first,AVG_num_of_mm_first_guide, AVG_of_mm])

def parse_all(path, outfile):
	f2 = open(outfile, 'w')
	for dir in os.listdir(path):
		if "HOM" in dir:
			f2.write(str(parse_file(os.path.join(path,os.path.join(dir, "output.txt")))[1:-1])+ '\n')
	f2.close()

if __name__ == "__main__":
	parse_all("groups/itay_mayrose/shiranabad/CRISPR/gal/gals_alg/")