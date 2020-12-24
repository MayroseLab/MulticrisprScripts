import pickle
import Candidate
import subgroup_res


def sub_tree_display(path, candidates_lst, f, consider_homology, counter):
	#candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))

	def chagne_mm_to_lowercase(target_str, mm_lst):
		'''mm is mismatch'''
		target_in_lst = list(target_str)
		for place in mm_lst:
			target_in_lst[place] = target_in_lst[place].lower()
		return ''.join(target_in_lst)

	#filepath = path + "/the_table.html"
	#f = open(filepath, 'w')
	header_row = "<tr>\n\t<th></th>\n\t<th>sgRNA</th>\n\t<th>Sum of genes score</th>\n\t<th>Genes</th>\n\t<th>Genes gcore</th>\n\t<th>Target site</th>\n</tr>\n"
	if consider_homology == True:
		f.write("<table id="+str(counter)+" style='display:none; width:100%'>\n")
	else:
		f.write("<table>\n")  	
		f.write(header_row)
	
	for c in candidates_lst:
		f.write("<tr>\n")
	#    print("candidate:\n",c)
		num_of_targets = 0
		for targets in c.targets_dict.values():
			num_of_targets += len(targets)
	#    print("num_of_targets: ",num_of_targets)
		first_target = 1
		first_gene = 1
		l = list(c.targets_dict.items())
		l.sort(key = lambda item: c.genes_score_dict[item[0]], reverse = True)
		#print(l)
		#for gene, targets in c.targets_dict.items():
		for gene, targets in l:
			#print("gene: ", gene)
			#print("targets: ",targets)

			first_target = 1
			for target in targets:
				if first_target == 1 and first_gene==1:
					f.write("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"' onclick='toggle(\""+c.seq+"\")'>+</td>\n")
					f.write("\t<td class='"+c.seq+"' onclick='toggle(\""+c.seq+"\")'>+</td>\n")
					f.write("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"'>"+c.seq+"</td>\n")
					f.write("\t<td class='"+c.seq+"'>"+c.seq+"</td>\n")
					f.write("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"'>"+str(c.cut_expectation)[:5]+"</td>\n")
					f.write("\t<td class='"+c.seq+"'>"+str(c.cut_expectation)[:5]+"</td>\n")
					f.write("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(len(targets))+"'>"+gene+"</td>\n")
					f.write("\t<td class='"+c.seq+"'>"+gene+"</td>\n")
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]
					#for Gene, score in c.genes_score_dict.items():
					#	if Gene == gene:
					f.write("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(len(targets))+"'>"+score+"</td>\n")
					f.write("\t<td class='"+c.seq+"'>"+score+"</td>\n")
							#break
					f.write("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					f.write("</tr>\n")
					first_target = 0
					continue
				if first_target != 1:
					f.write("<tr class='"+c.seq+"_r' style='display:none'>\n")
					f.write("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					f.write("</tr>\n")
				if first_target == 1 and first_gene != 1:
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]
					f.write("<tr class='"+c.seq+"_r' style='display:none'>\n")
					f.write("\t<td rowspan='"+str(len(targets))+"'>"+gene+"</td>\n")
					f.write("\t<td rowspan='"+str(len(targets))+"'>"+score+"</td>\n")
					f.write("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					f.write("</tr>\n")
					first_target = 0
			first_gene = 0
			



	f.write("</table>\n")
	#f.write("<script>\nfunction toggle(sgRNA) {\n\t//document.getElementById(sgRNA).style.display = document.getElementById(sgRNA).style.display === 'none' ? 'table-row' : 'none';\n    var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	#f.write("        lst_r[i].style.display = lst_r[i].style.display === 'none' ? 'table-row' : 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	#f.write("    }\n  }\n </script>")


	#f.close()
def genes_of_sub_candidates_lst(sub_candidates_lst):
	res = set()
	for c in sub_candidates_lst:
		for gene in c.genes_score_dict.keys():
			res.add(gene)
	return list(res)

def tree_display(path, consider_homology = False):
	candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))
	filepath = path + "/the_table.html"
	f = open(filepath, 'w')
	f.write("The designed sgRNAs for the genes in your input are listed in the table below.<br>Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")
	f.write("<table style=\"width:70%; font-size:20px; font-family:aharoni;\">\n")
	f.write("  <tr>\n    <th><b>Homologous subgroup:</b> Genes</th>\n  </tr>\n")
	f.write("</table>\n")
	counter = 0;
	if consider_homology == False:
		sub_tree_display(path, candidates_lst, f, consider_homology,counter)
	else:
		#for sub_candidates_lst in candidates_lst:
		for subgroup_item in candidates_lst:

			# create the main table
			counter +=1
			#genes = str(genes_of_sub_candidates_lst(sub_candidates_lst))[1:-1]
			genes = str(subgroup_item.genes_lst)[1:-1].replace("'","")

			f.write("<table style=\"width:70%; font-size:20px; font-family:aharoni; \">\n")
			f.write("  <tr>\n       <td onclick='show("+str(counter)+")'>+ "+subgroup_item.name+ ": "+genes+"</td>\n   </tr>\n") #why the + before "+genes+" ?
			f.write("</table>\n")
			sub_tree_display(path, subgroup_item.candidate_lst, f, consider_homology, counter)
		counter = 0;	
		
	f.write("<script>\nfunction toggle(sgRNA) {\n\t//document.getElementById(sgRNA).style.display = document.getElementById(sgRNA).style.display === 'none' ? 'table-row' : 'none';\n    var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	f.write("        lst_r[i].style.display = lst_r[i].style.display === 'none' ? 'table-row' : 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	f.write("    }\n  }\n ")
	f.write("function show(counter) {\n document.getElementById(counter).style.display = document.getElementById(counter).style.display === 'none' ? 'block' : 'none';\n  }\n")
	f.write("</script>\n")
	f.close()
#test

if __name__ == "__main__":
	print("ASDF")
	tree_display("/groups/itay_mayrose/galhyams/EData/red_sol_fam/output_homology")
