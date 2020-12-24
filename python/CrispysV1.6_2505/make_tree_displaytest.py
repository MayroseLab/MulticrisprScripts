import pickle
import Candidate


def tree_display(path):
	candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))

	def chagne_mm_to_lowercase(target_str, mm_lst):
		'''mm is mismatch'''
		target_in_lst = list(target_str)
		for place in mm_lst:
			target_in_lst[place] = target_in_lst[place].lower()
		return ''.join(target_in_lst)

	filepath = path + "\\the_table.html"
	f = open(filepath, 'w')
	header_row = "<tr>\n\t<th></th>\n\t<th>sgRNA</th>\n\t<th>scoret</th>\n\t<th>genes</th>\n\t<th>Genes score</th>\n\t<th>Target site</th>\n</tr>\n"
	print("<table>\n")
	print(header_row)
	for c in candidates_lst:
		print("<tr>\n")
	#    print("candidate:\n",c)
		num_of_targets = 0
		for targets in c.targets_dict.values():
			num_of_targets += len(targets)
	#    print("num_of_targets: ",num_of_targets)
		first_target = 1
		first_gene = 1
		for gene, targets in c.targets_dict.items():
			first_target = 1
			for target in targets:
				if first_target == 1 and first_gene==1:
					print("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"' onclick='toggle(\""+c.seq+"\")'>+</td>\n")
					print("\t<td class='"+c.seq+"' onclick='toggle(\""+c.seq+"\")'>+</td>\n")
					print("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"'>"+c.seq+"</td>\n")
					print("\t<td class='"+c.seq+"'>"+c.seq+"</td>\n")
					print("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(num_of_targets)+"'>"+str(c.cut_expectation)+"</td>\n")
					print("\t<td class='"+c.seq+"'>"+str(c.cut_expectation)+"</td>\n")
					print("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(len(targets))+"'>"+gene+"</td>\n")
					print("\t<td class='"+c.seq+"'>"+gene+"</td>\n")
					#score = c.genes_score_dict[gene]
					for Gene, score in c.genes_score_dict.items():
						if Gene == gene:
					print("\t<td class='"+c.seq+"' style='display:none' rowspan='"+str(len(targets))+"'>"+str(score)+"</td>\n")
					print("\t<td class='"+c.seq+"'>"+str(score)+"</td>\n")
							break
					print("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					print("</tr>\n")
					first_target = 0
					continue
				if first_target != 1:
					print("<tr class='"+c.seq+"_r' style='display:none'>\n")
					print("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					print("</tr>\n")
				if first_target == 1 and first_gene != 1:
					#score = c.genes_score_dict[gene]
					print("<tr class='"+c.seq+"_r' style='display:none'>\n")
					print("\t<td rowspan='"+str(len(targets))+"'>"+gene+"</td>\n")
					print("\t<td rowspan='"+str(len(targets))+"'>"+str(score)+"</td>\n")
					print("\t<td>"+chagne_mm_to_lowercase(target[0], target[1].keys())+"</td>\n")
					print("</tr>\n")
					first_target = 0
			first_gene = 0
			



	print("</table>\n")
	print("<script>\nfunction toggle(sgRNA) {\n\t//document.getElementById(sgRNA).style.display = document.getElementById(sgRNA).style.display === 'none' ? 'table-row' : 'none';\n    var lst_r = document.getElementsByClassName(sgRNA.concat('_r'));\n    for(var i = 0; i < lst_r.length; ++i) {\n")
	print("        lst_r[i].style.display = lst_r[i].style.display === 'none' ? 'table-row' : 'none';\n    }\n    var lst = document.getElementsByClassName(sgRNA);\n    for(var i = 0; i < lst.length; ++i) {\n        lst[i].style.display = lst[i].style.display === 'none' ? 'table-cell' : 'none';")
	print("    }\n  }\n </script>")


	f.close()
	
#test

if __name__ == "__main__":
	tree_display("/groups/itay_mayrose/galhyams/EData/AT5G14720/E")
	print("ASDFASDFASFASDF")
