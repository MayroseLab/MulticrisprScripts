import pickle
import Candidate
import subgroup_res
import CasSites


def sub_tree_display(path, candidates_lst, f, consider_homology, counter, genes_names, genes_lst):
	#candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))

	def chagne_mm_to_lowercase(target_str, mm_lst):
		'''mm is mismatch'''
		target_in_lst = list(target_str)
		for place in mm_lst:
			target_in_lst[place] = target_in_lst[place].lower()
		return ''.join(target_in_lst)

	#filepath = path + "/the_table.html"
	#f = open(filepath, 'w' 

	header_row = "sgRNA index,sgRNA,Score,Genes,Genes score,Target site,#mms,Position\n"#new
	#header_row = "<tr>\n,<th></th>\n,<th>sgRNA</th>\n,<th>Score</th>\n,<th>Genes</th>\n,<th>Genes score</th>\n,<th>Target site</th>\n,<th>#mms</th>\n,<th>Find off-targets</th>\n</tr>\n"#new

	if consider_homology == True:
		f.write("table id="+str(counter)+"\n")#new
	#else: 
		#f.write("Filter out targets with more than\n<select id=\"selectMM\" onchange=\"javascript:filterMM();\">\n")#new
		#f.write(",<option value=\""+str(100)+"\">""</option>\n")#new
		#for i in range(1,21):#new
		#	f.write(",<option value=\""+str(i)+"\">"+str(i)+"</option>\n")#new
		#f.write("</select>\n mismatches\n")#new
		#f.write("<table class=\"table\">\n")  	#new
	f.write(header_row)
	sgRNA_index = 0
	for c in candidates_lst:
		sgRNA_index += 1
		#f.write("<tr class='top_row'>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n</tr>\n")
		#f.write("<tr>\n")      
	#    print("candidate:\n",c)
		num_of_targets = 0
		for targets in c.targets_dict.values():
			num_of_targets += len(targets)
	#    print("num_of_targets: ",num_of_targets)
		first_target = 1
		first_gene = 1
		l = list(c.targets_dict.items())
		l.sort(key = lambda item: c.genes_score_dict[item[0]], reverse = True)
		#for gene, targets in c.targets_dict.items():
		for gene, targets in l:
			targets.sort(key = lambda target: len(target[1]))#Galll!!! sort the targets by number of mms 
			gene_seq = genes_lst[genes_names.index(gene)]
			seen_sites = dict()
			first_target = 1
			for target in targets:
				pos = find_pos(target,gene_seq, seen_sites)
				#print ('pos:', pos)
				if first_target == 1 and first_gene==1:
					#f.write(","+c.seq)#+\n")
					#f.write(", class='"+c.seq+"' onclick='show_row(\""+c.seq+"\")'>+</td>\n")
					#f.write(", class='"+c.seq+"' style='display:none'>"+c.seq+"</td>\n")
					#f.write(", class='"+c.seq+"'>"+c.seq+"</td>\n")
					f.write(str(sgRNA_index) + '.,' + c.seq+","+str(c.cut_expectation)[:5])#+"\n")
					#f.write(", class='"+c.seq+"'>"+str(c.cut_expectation)[:5]+"</td>\n")
					f.write(","+gene)#+"</td>\n")
					#f.write("," class='"+c.seq+"'>"+gene+"</td>\n")
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]
					#for Gene, score in c.genes_score_dict.items():
					#	if Gene == gene:
					f.write(","+score)#+"\n")
					#f.write(",<td class='"+c.seq+"'>"+score+"</td>\n")
					f.write(","+chagne_mm_to_lowercase(target[0], target[1].keys()))#+"\n")
					#print ("target[0] " , target[0] , " target[1].keys() " , target[1] ,"len: ", len(target[1]),"\n")#new
					f.write(","+str(len(target[1])))#+"\n")#new
					f.write(","+ pos)# +"\n")#new
					#f.write(",<td class='"+c.seq+"' style='display:none'><a target=\"_blank\" href=\"http://tefor.net/crispor/crispor.cgi?seq="+c.seq+"&pam=NGG\">CRISPOR</a> | <a target=\"_blank\" href=\"http://crista.tau.ac.il/findofftargets.html#form/sgRNA_seq="+c.seq+"\">CRISTA</a></td>\n")
					#f.write(",<td class='"+c.seq+"'><a target=\"_blank\" href=\"http://tefor.net/crispor/crispor.cgi?&seq="+c.seq+"&pam=NGG\">CRISPOR</a> | <a target=\"_blank\" href=\"http://crista.tau.ac.il/findofftargets.html#form/sgRNA_seq="+c.seq+"\">CRISTA</a></td>\n")
					f.write("\n")
					first_target = 0
					continue
				if first_target != 1:
					#f.write("<tr class='"+c.seq+"_r' style='display:none'>\n")
					f.write(str(sgRNA_index) + ".,,,,,")
					#f.write(",<td></td>\n")
					#f.write(",<td></td>\n")
					#f.write(",<td></td>\n")
					#f.write(",<td></td>\n")
					f.write(chagne_mm_to_lowercase(target[0], target[1].keys()))#+"</td>\n")
					f.write(","+str(len(target[1])))#+"\n")#new
					f.write(","+pos)#+"</td>\n")#new
					#f.write(",<td></td>\n")
					f.write("\n")
				if first_target == 1 and first_gene != 1:
					f.write(str(sgRNA_index) + ".,,,")
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]
					#f.write("<tr class='"+c.seq+"_r' style='display:none'>\n")
					#f.write(",<td></td>\n")
					#f.write(",<td></td>\n")
					#f.write(",<td></td>\n")
					f.write(gene)#+"</td>\n")
					f.write(","+score)#+"</td>\n")
					f.write(","+chagne_mm_to_lowercase(target[0], target[1].keys()))#+"</td>\n")
					f.write(","+str(len(target[1])))#+"</td>\n")#new
					f.write(","+pos)#+"</td>\n")#new
					f.write("\n")
					#f.write("</tr>\n")
					first_target = 0
			first_gene = 0
		#f.write("<tr class='top_row'>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n</tr>")	


	#f.write("<tr class='top_row'>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n,<td></td>\n</tr>")
	#f.write("</table>\n")
	#f.write("<script>\nfunction toggle(sgRNA) {\n,//document.getElementById(sgRNA).style.display = document.getElementById(sgRNA).style.display === 'none' ? 'table-row' : 'none';\n    var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	#f.write("        lst_r[i].style.display = lst_r[i].style.display === 'none' ? 'table-row' : 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	#f.write("    }\n  }\n </script>")
	#f.close()
	
def genes_of_sub_candidates_lst(sub_candidates_lst):
	res = set()
	for c in sub_candidates_lst:
		for gene in c.genes_score_dict.keys():
			res.add(gene)
	return list(res)

def test_find_pos():
	target, gene_sequence = ['AAAAAAAAAAAAAAAAAAGG'], 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGG'
	seen_sites = dict()
	position = 0
	for i in range(99):
		#r = find_pos(target, gene_sequence, seen_sites)
		position = gene_sequence.find(target[0], position+1)
		if position == -1:
			break
	r = find_pos(target, gene_sequence, seen_sites)


	
def find_pos(target, gene_sequence, seen_sites):
	'''sgRNA_targets is a list of target sites
	returns'''
	#seen_sites = dict()
	#res = dict() #key:targets val: position in gene
	#for target in sgRNA_targets:
	target_seq = target[0]
	if target_seq in seen_sites:
		directions_lst = seen_sites[target_seq]
	else:
		directions_lst = [0,0]
	position = gene_sequence.find(target_seq, directions_lst[0])#+ '+' #item[0] is the target site seq
	if position != -1:
		update_seen_sites_dict(seen_sites, target_seq, 0, position)
		position = str(position) + '+'
	#update_positions_dict(seen_sites, target ,sg_position)
	else:
		position = CasSites.give_complementary(gene_sequence).find(target_seq, directions_lst[1])# + "-"
		update_seen_sites_dict(seen_sites, target_seq, 1, position)
		position = str(position) + '-'
	if position == -1:
		position = ''
		#update_positions_dict(res, target ,sg_position)
		#update_seen_sites_dict(seen_sites, target_seq, 1)
	return position

def update_positions_dict(d, site_seq, position):
	if site_seq in d:
		d[site_seq] = d[site_seq] + [position]
	else:
		d[site_seq] = [position]


def update_seen_sites_dict(d, site_seq, direction, position):
	'''
	d: dict: key: site_seq; val: [num of occurrences, directions_lst]
	direction: 0 for sense, 1 for antisence
	'''
	if site_seq in d:
		directions_lst = d[site_seq]
		directions_lst[direction] = position + 1
		d[site_seq] = directions_lst
	else:
		directions_lst = [0,0]
		directions_lst[direction] = position + 1
		d[site_seq] = directions_lst

		
def tree_display(path, consider_homology = False, set_cover = False):
	if set_cover:
		candidates_lst = pickle.load(open(path + "/greedy_cover.p", "rb"))		
	else:
		candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))
	genes_names = pickle.load(open(path + "/genesNames.p", "rb"))
	genes_list = pickle.load(open(path + '/genesList.p', 'rb'))

	filepath = path + "/CRISPys_output.csv"
	f = open(filepath, 'w')
	counter = 0;
	if consider_homology == False:
		sub_tree_display(path, candidates_lst, f, consider_homology,counter, genes_names, genes_list)
	else:
		f.write("The designed sgRNAs for the genes in your input are listed in the table below. Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")
		#f.write("Filter out targets with more than\n<select id=\"selectMM\" onchange=\"javascript:filterMM();\">\n")#new
		#f.write(",<option value=\""+str(100)+"\">""</option>\n")#new
		#for i in range(1,21):#new
		#	f.write(",<option value=\""+str(i)+"\">"+str(i)+"</option>\n")#new
		#f.write("</select>\n mismatches\n")#new
		#f.write("<table style=\"font-size:20px; font-family:aharoni;\">\n")
		#f.write("  <tr>\n    <th><b>Homologous subgroup:</b> Genes</th>\n  </tr>\n")
		#f.write("</table>\n")
		#for sub_candidates_lst in candidates_lst:
		for subgroup_item in candidates_lst:

			# create the main table
			counter +=1
			#genes = str(genes_of_sub_candidates_lst(sub_candidates_lst))[1:-1]
			genes = str(subgroup_item.genes_lst)[1:-1].replace("'","")

			#f.write("<table class=\"table group\" style=\"font-size:20px; font-family:aharoni; \">\n")
			#f.write("  <tr class=\"group\">\n       <td onclick='show("+str(counter)+")'>+ "+subgroup_item.name+ ": "+genes+"</td>\n       <td style='display:none' onclick='show("+str(counter)+")'>- "+subgroup_item.name+ ": "+genes+"</td>\n   </tr>\n") #why the + before "+genes+" ?
			#f.write("</table>\n")
			sub_tree_display(path, subgroup_item.candidate_lst, f, consider_homology, counter, genes_names, genes_list)
		counter = 0;	
		
	#f.write("<script>\nfunction hide(sgRNA) {\n, var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	#f.write("        if (lst_r[i].classList.contains('hideFromFilter')){lst_r[i].classList.remove('hideFromFilter');}\n        lst_r[i].style.display = 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	#f.write("    }\n  }\n ")
	#f.write("function show_row(sgRNA) {\n, var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n , var mm = document.getElementById(\"selectMM\").value;    for(var i = 0; i < lst_r.length; ++i) {\n")
	#f.write("        var myMM = lst_r[i].children[\"mms\"].innerText;\n ,if (parseInt(myMM)  <= parseInt(mm)) {\n lst_r[i].style.display = 'table-row';\n,}\n ,else {lst_r[i].classList.add('hideFromFilter');}\n   }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	#f.write("    }\n  }\n ")
	#f.write("function show(counter) {\n document.getElementById(counter).style.display = document.getElementById(counter).style.display === 'none' ? 'block' : 'none';\n  }\n")
	#f.write("function filterMM(){\n var mm = document.getElementById(\"selectMM\").value;\n var tables = document.getElementsByClassName(\"table\");\n for (var t = 0; t < tables.length; ++t)\n {\n ,table = tables[t];\n ,tr = table.getElementsByTagName(\"tr\");\n ,for (i = 0; i < tr.length; i++)\n ,{\n    if(tr[i].style.display === 'none'){\n      if (tr[i].classList.contains('hideFromFilter')){}\nelse{continue;}\n    } ,,var row = tr[i].children[\"mms\"];\n ,, if (row){\n ,,, var myMM = row.innerText;\n ,,, if (parseInt(myMM)  > parseInt(mm)) {\n        tr[i].classList.add('hideFromFilter');\n        tr[i].style.display = \"none\";\n       }\n ,,, else{tr[i].style.display = \"table-row\";tr[i].classList.remove('hideFromFilter');}\n ,,}\n ,}\n }}\n")
	#f.write("</script>\n")
	f.close()
#test

if __name__ == "__main__":
#	print("ASDF")
	#test_find_pos()
	tree_display("/groups/itay_mayrose/galhyams/1516893877")
