from __future__ import print_function

import sys
import numpy as np
import Candidate

import cplex
from cplex.exceptions import CplexError

import math



def populatebyrow(prob, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows):
	prob.objective.set_sense(prob.objective.sense.minimize)

	prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype,
					   names=my_colnames)



	prob.linear_constraints.add(lin_expr=rows, senses=my_sense, rhs=my_rhs, names=my_rownames)

def mipex1(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows):

	try:
		my_prob = cplex.Cplex()

		if pop_method == "r":
			handle = populatebyrow(my_prob, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)
		elif pop_method == "c":
			handle = populatebycolumn(my_prob)
		elif pop_method == "n":
			handle = populatebynonzero(my_prob)
		else:
			raise ValueError('pop_method must be one of "r", "c" or "n"')

		my_prob.solve()
	except CplexError as exc:
		print(exc)
		return exc

	#print()
	# solution.get_status() returns an integer code
	#print("Solution status = ", my_prob.solution.get_status(), ":", end=' ')
	# the following line prints the corresponding string
	#print(my_prob.solution.status[my_prob.solution.get_status()])
	#print("Solution value  = ", my_prob.solution.get_objective_value())

	numcols = my_prob.variables.get_num()
	numrows = my_prob.linear_constraints.get_num()

	slack = my_prob.solution.get_linear_slacks()
	x = my_prob.solution.get_values()

	#for j in range(numrows):
	#	print("Row %d:  Slack = %10f" % (j, slack[j]))
	#for j in range(numcols):
	#	print("Column %d:  Value = %10f" % (j, x[j]))
	return my_prob
		

def CplexSetCover(candidatesLst, genesLst, thr):
	eps = 0.001
	my_obj = [1.0 for i in range(len(candidatesLst))]
	my_ub = [1.0 for i in range(len(candidatesLst))]
	my_lb = [0.0 for i in range(len(candidatesLst))]
	my_ctype = "I" * len(candidatesLst)
	my_colnames = ["x" + str(i+1) for i in range(len(candidatesLst))]

	my_rhs = [1.0 for i in range(len(genesLst))]
	my_rownames = ["r" + str(i+1) for i in range(len(genesLst))]
	my_sense = "G" * len(genesLst) #upper bound

	rows = make_rows_set_cover(candidatesLst, genesLst, thr - eps)

	return mipex1("r", my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)#(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)

def Cplex_fuzzy_set_cover(candidatesLst, genesLst, thr, FULL_COVER_THR = 0.99):
	eps = 0.001
	# data common to all populateby functions
	#first, make a table
	#candidatesLst = candidatesLst[:500]#test
	#print(candidatesLst)
	#genesLst = genesLst[:5]
	#table_A = make_a_table(candidatesLst, genesLst, thr)
	#print(table_A)
	#second, define the variables
	thr = -math.log2(1-thr)#can be math.log as well
	my_obj = [1.0 for i in range(len(candidatesLst))]
	my_ub = [1.0 for i in range(len(candidatesLst))]
	my_lb = [0.0 for i in range(len(candidatesLst))]
	my_ctype = "I" * len(candidatesLst)
	my_colnames = ["x" + str(i+1) for i in range(len(candidatesLst))]
	#my_colnames = ["x1", "x2", "x3", "x4"]
	my_rhs = [thr for i in range(len(genesLst))]
	my_rownames = ["r" + str(i+1) for i in range(len(genesLst))]
	my_sense = "G" * len(genesLst) #change for upper bound somehow
	#input matrix for the problem
	rows = make_rows_fuzzy_set_cover(candidatesLst, genesLst, thr - eps, 1-eps)
	#rows = []
	#for i in range(len(genesLst)):
	#    rows.append([["x" + str(i+1)], [1.0]])
	#print(rows)
	#print("\n\n\n")
	#rows = [["x" + str(i)], [1.0] for i in range len(genesLst)]

	#rows = [[["x1", "x2", "x3", "x4"], [-1.0, 1.0, 1.0, 10.0]],
			#[["x1", "x2", "x3"], [1.0, -3.0, 1.0]],
			#[["x2", "x4"], [1.0, -3.5]]]

	#print(my_colnames)
	#print(rows)
	prob_obj = mipex1("r", my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)#(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)
	#print("prob obj\n\n\n\n\n\n", prob_obj)
	return prob_obj

	
def make_rows_set_cover(candidatesLst, genesLst, thr):
	'''make a table of a_ij: does the i gene is covered by the j sgRNA'''
	#A = np.zeros((len(genesLst), len(candidatesLst)), dtype=bool)
	rows = []
	for i in range(len(genesLst)): # a row
		row = [[],[]]#row = list(list())
		for j in range(len(candidatesLst)):
			if genesLst[i] in candidatesLst[j].genes_score_dict:
				if candidatesLst[j].genes_score_dict[genesLst[i]]> thr:
					row[0].append('x' + str(j+1))
					row[1].append(1.0)
					#print("gene " + str(i)+ "is covered by candidate" + str(j))
					#A[i][j] = True
		rows.append(row)
	return rows
	
def make_rows_fuzzy_set_cover(candidatesLst, genesLst, thr, FULL_COVER_THR):
	'''make a table of a_ij: does the i gene is covered by the j sgRNA'''
	#A = np.zeros((len(genesLst), len(candidatesLst)), dtype=bool)
	rows = []
	for i in range(len(genesLst)): # a row
		row = [[],[]]#row = list(list())
		for j in range(len(candidatesLst)):
			if genesLst[i] in candidatesLst[j].genes_score_dict:
				#if candidatesLst[j].genes_score_dict[genesLst[i]]> thr:
				row[0].append('x' + str(j+1))
				#print("ASDFASDFASDFASDFD")
				#print(candidatesLst[j].genes_score_dict[genesLst[i]])
				if candidatesLst[j].genes_score_dict[genesLst[i]] == 1:
					row[1].append(-math.log2(1-(FULL_COVER_THR)))
				else:
					#print("score: ", candidatesLst[j].genes_score_dict[genesLst[i]])
					row[1].append(-math.log2(1-(candidatesLst[j].genes_score_dict[genesLst[i]])))
				#print("gene " + str(i)+ "is covered by candidate" + str(j))
					#A[i][j] = True
		rows.append(row)
	return rows	
	
def make_a_table(candidatesLst, genesLst, thr):
	'''make a table of a_ij: does the i gene is covered by the j sgRNA'''
	A = np.zeros((len(genesLst), len(candidatesLst)), dtype=bool)

	for j in range(len(candidatesLst)):
		for i in range(len(genesLst)):
			if genesLst[i] in candidatesLst[j].genes_score_dict:
				if candidatesLst[j].genes_score_dict[genesLst[i]]> thr:
					A[i][j] = True
	return A

