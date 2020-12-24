__author__ = 'ItayM5'
import re ##used to be regex
import regex

def get_sites(gene, min_length=20, max_length=20, start_with_G = False, where_in_gene = 1, CRISTA = False, PAMs = ["GG"]):
	'''
	:param gene:
	:param min_length:
	:param max_length:
	:param start_with_G:
	:param where_in_gene: forword to this position the sgRNA are ignored
	:return:
	'''
	res = []
	if len(gene) < max_length+3:
		return res
	for length in range(min_length, max_length +1):
		for i in range(len(PAMs)):
			if (start_with_G):
				SiteAndPAM = "G" + "."*length + PAMs[i] #it is acually NGG
			else:
				SiteAndPAM = "."*(length +1) + PAMs[i] #it is acually NGG
			if CRISTA:
				SiteAndPAM = "."*3 + SiteAndPAM + '.'*3 #make sure that this is the needed input
			compiled = regex.compile(SiteAndPAM)
			where_in_gene = int(len(gene)*where_in_gene)
			founds_sense = regex.findall(compiled, gene[:where_in_gene], overlapped=True)
			founds_antisense = regex.findall(compiled, give_complementary(gene)[:where_in_gene], overlapped=True)
			if CRISTA:
				founds = founds_sense + founds_antisense
			else:
				founds = [seq[:-3] for seq in founds_sense if 'N' not in seq[:-3]] + [seq[:-3] for seq in founds_antisense if 'N' not in seq[:-3]]
			res += founds
	return res

def get_sites_test(gene, min_length=20, max_length=20, start_with_G=False, where_in_gene = 1):
	res = []
	SiteAndPAM = "."*(20 +1) + "GG" #it is acually NGG
	compiled = regex.compile(SiteAndPAM)
	#where_in_gene = int(len(gene)*where_in_gene)
	founds_sense = regex.findall(compiled, gene, overlapped=True)
	#print("gene", gene)
	#print("found sense", founds_sense)
	founds_antisense = regex.findall(compiled, give_complementary(gene), overlapped=True)
	founds = [seq[:-3] for seq in founds_sense] + [seq[:-3] for seq in founds_antisense]
	res = founds
	return res

def give_complementary(seq):
	res = []
	for i in range(len(seq)) :
		if seq[len(seq)-1-i] == 'A':
			res.append('T')
		elif seq[len(seq)-1-i] == 'T':
			res.append('A')
		elif seq[len(seq)-1-i] == 'C':
			res.append('G')
		elif seq[len(seq)-1-i] == 'G':
			res.append('C')
		elif seq[len(seq)-1-i] == 'N':
			res.append('N')
	return ''.join(res)

def give_complementary_old(seq):
	res = []
	for letter in seq:
		if letter == 'A':
			res.append('T')
		elif letter == 'T':
			res.append('A')
		elif letter == 'C':
			res.append('G')
		elif letter == 'G':
			res.append('C')
		elif letter == 'N':
			res.append('N')
	return ''.join(res)


def get_targets_sites_from_exons_lst(exons_lst, original_range_in_gene = [0,1], min_length= 20, max_length = 20,start_with_G = False, CRISTA = False, PAMs = ["GG"]):
	if original_range_in_gene[1] <= original_range_in_gene[0]:
		print("The range of the targts on the gene is not in the right format")
		exit(-1)
	if max_length < min_length:
		print("The range of the lengths of the sgRNA is not in the right format")
		exit(-1)
	res = []
	lengths = list(map(lambda x: len(x), exons_lst))
	gene_length = sum(lengths)
	range_in_gene = list(map(lambda x: int(x * gene_length), original_range_in_gene))
	exons_lst = list(map(lambda seq: seq.upper(), exons_lst)) #converting to upper-case
	where_in_gene = 1
	for i in range(1, len(lengths)):
		lengths[i] = lengths[i-1] + lengths[i]
	for i in range(len(exons_lst)):
		if i == 0:
			if range_in_gene[0] < lengths[i]:
				#if range_in_gene[1]*gene_length > lengths[i]:
				res += get_sites(exons_lst[i][range_in_gene[0] : min(lengths[i], range_in_gene[1])], min_length, max_length, start_with_G, where_in_gene, CRISTA, PAMs)
		elif max(range_in_gene[0], lengths[i-1]) < min(lengths[i], range_in_gene[1]):
			res += get_sites(exons_lst[i][max(range_in_gene[0]  - lengths[i-1], 0) : min(lengths[i] - lengths[i-1], range_in_gene[1] - lengths[i-1])], min_length, max_length, start_with_G, where_in_gene, CRISTA, PAMs)
	return res

def test_2():
	gene = ["TTTATGTCAACTTTTTCAATCTAATAGATCAATGAATTGTAAACTTTTTTCGACCACAAAATGATGCTTCCAAATACAAACAAAACCTGATGCAATCAGTCAATACCTTCCAACTTTAGAACACATATATGTAGCAATGCTCCTACAGTTTACTTTTCTATCTTTTAGCCTAATCATTTACTCTCATATTTTTTCTTTAAACTAGAAAGTTCAGAATCCAAATATAATATCATCTCCTTCTCTCTATTACAGCAATGGTTTTGGTTGATAACCATGCTGGAAAAGATGGTGCAGAAGATGGTAATATGGTTGATTTTCGAGGAAATCCGGTGGATAAGTCTAGGACAGGGGGATGGCTAGCTGCAGGACTTATCCTAGGAACTGAGCTATCAGAAAGGGTATGTGTTATGGGGATTTCGATGAATTTAGTGACGTACTTAGTTGGAGATTTACATCTTCCATCCTCCAAATCTGCCAACATTGTCACCAATTTCATGGGGACACTTAATCTTCTTGGTCTTCTAGGTGGTTTCTTGGCAGATGCTAAACTCGGACGTTATCTGACTGTTGGAATCTTTGCTTCAATTGCTGCTGTGGGGGTTACGCTTTTGACATTGGCGACATCCATTCCAGGCATGAAGCCGCCTGAATGTAACCCAAGAAAAAGTGGTCACTGCATTGAAGCCAGTGGCCAGCAGCTTGCTCTTCTCTATACGGCGCTTTACATCCTAGCTCTTGGTGGTGGTGGAATTAAGTCAAATGTCTCCGGGTTTGGTTCAGACCAATTTGACTCATCAGATCCTAAGGAGAACAAGTCCATGATATACTTCTTCAACAGATTCTATTTCTGCATAAGCCTTGGTTCTCTGTTTGCAGTGACTGTGCTGGTGTACTTACAAGACAATGTAGGAAGAGGATGGGGATATGGGATATCAGCAGGCACAATGGTCCTCGGGGTCGCTGTATTGATTGGTGGAACGACGTTGTATCGATTCAAGAAGCCTCAAGGAAGTCCTTTGACTATCATATGGAGGGTTCTGCTTTTAGCTTGGAGGAAGAGAAAGCTTAGTTACCCTTCTGATACTGGCTTCTTGAATGAATATCACAATGCCAAAGTCCCACATACACATATGTTGAGGTGTCTTGACAAGGCAGCCATTCTTGATGACTCTGCAGCTGCAAATGAGAATAGCAAGAATCGTTGGATAGTTTCAACAGTTACAGAAGTCGAAGAAGTGAAAATGGTGCTCAAATTGATTCCCATATGGTCCACATGCATACTTTTTTGGACAGTATACTCTCAGATGAATACCTTCACCATTGAACAAGCTACCTTCATGAACCGGAATGTTGGAAACTTTGCTGTCCCTGCAGGTTCCTTATCCGTGTTTCTCTTTATTAGCATACTTCTGTTTACTTCCATAAACGAAAGGGTCACAGTTCGTATTGCCAGAAAAATCACTCACAACAGCCAAGGAATCACAAGCCTTCAGAGAGTTGGAATTGGACTACTACTCTCTATTGTTGGTATGGTAGCTTCAGCTCTGGTAGAAAAACGACGAAGGGAACATGCCATCCATCATAACTTCAAGATAAGCGCGTTTTGGTTAGTGCCTCAATTCTTCATTGTAGGTGCTGGGGAAGCTTTTGCCTATGTAGGACAGCTAGAGTTTTTCATCAGGGAGGCACCAGAAGGGATGAAATCTATGAGCACAGGCCTATTTCTCAGCACACTCTCGATGGGATATTTCGTGAGTAGTTTGCTAGTATTCGTTGTACAGAAAGCAACAAAAGGAAGATGGCTTAAAAGCAATTTAAACAAAGGAAAACTGGATTTATTCTACTGGTTGCTAGCAGTTCTCGGAGTAATTAATTTCTTGATTTTCATTGCATTTTCAATGAAACACCAATACAAGGTGCAGAAACTTAGCAGTATTGAGGATTCTGCAGAAGAGCTCGGGAGTTGGAAGGATTTGACCCTCGACAACAAGGAAAAGAAACTCGAAGCAGACGAGAAGGTGGAAGCTTAAATACAGCATATTAGCTTTCAATGAATCATTCATTTCCAGAGTTTGTAATATAGAACCGTATTCAATTATCAAAGACGTCAATACAAATTTGCTACCAGTCTTGAGTTCTGTTTAGATTAAAACCTTGGATATTAGAGTGCAGAAATATGATCAATTCAGAAAGATATTTACACTTCAAATTCTCACTAAA"]
	print(get_targets_sites_from_exons_lst(gene))
	print(get_sites("".join(gene)))

if __name__ == "__main__":
	test_2()