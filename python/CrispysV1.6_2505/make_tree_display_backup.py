import pickle
import Candidate


def tree_display(path, homology = False):
	candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))

	def chagne_mm_to_lowercase(target_str, mm_lst):
		'''mm is mismatch'''
		target_in_lst = list(target_str)
		for place in mm_lst:
			target_in_lst[place] = target_in_lst[place].lower()
		return ''.join(target_in_lst)

	filepath = path + "/the_table.html"
	f = open(filepath, 'w')
	header_row = "<tr>\n\t<th></th>\n\t<th>sgRNA</th>\n\t<th>Sum of genes score</th>\n\t<th>Genes</th>\n\t<th>Genes gcore</th>\n\t<th>Target site</th>\n</tr>\n"
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
	f.write("<script>\nfunction toggle(sgRNA) {\n\t//document.getElementById(sgRNA).style.display = document.getElementById(sgRNA).style.display === 'none' ? 'table-row' : 'none';\n    var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	f.write("        lst_r[i].style.display = lst_r[i].style.display === 'none' ? 'table-row' : 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	f.write("    }\n  }\n </script>")


	f.close()


	
#test

if __name__ == "__main__":
	print("ASDF")
	tree_display("/groups/itay_mayrose/galhyams/EData/red_sol_fam/output_homology")
