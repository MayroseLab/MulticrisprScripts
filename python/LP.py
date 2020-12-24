from __future__ import print_function
from Candidate import *
from Metric import cfd_funct
import numpy as np
from functools import reduce

def test_Cplex():
    import sys

    import cplex
    from cplex.exceptions import CplexError

    # data common to all populateby functions
    my_obj = [1.0, 2.0, 3.0, 1.0]
    my_ub = [40.0, cplex.infinity, cplex.infinity, 3.0]
    my_lb = [0.0, 0.0, 0.0, 2.0]
    my_ctype = "IIII"
    my_colnames = ["x1", "x2", "x3", "x4"]
    my_rhs = [20.0, 30.0, 0.0]
    my_rownames = ["r1", "r2", "r3"]
    my_sense = "LLE"
    prob = cplex.Cplex()

    prob.objective.set_sense(prob.objective.sense.maximize)

    prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype,
                       names=my_colnames)

    rows = [[["x1", "x2", "x3", "x4"], [-1.0, 1.0, 1.0, 10.0]],
            [["x1", "x2", "x3"], [1.0, -3.0, 1.0]],
            [["x2", "x4"], [1.0, -3.5]]]

    prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)

    prob.solve()



def test_Cplex1():
    '''
    :return: running a simple {0,1} linear program
    '''
    import sys
    import cplex
    from cplex.exceptions import CplexError

    prob = cplex.Cplex()
    # data common to all populateby functions
    my_obj = [1.0, 2.0, 3.0, 1.0]
    my_ub = [1.0, 1.0, 1.0, 1.0]
    my_lb = [0.0, 0.0, 0.0, 0.0]
    my_ctype = "IIII"
    my_colnames = ["x1", "x2", "x3", "x4"]
    my_rhs = [2.0]
    my_rownames = ["r1"]
    my_sense = "LLE"

    prob.objective.set_sense(prob.objective.sense.maximize)

    prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype,
                       names=my_colnames)

    rows = [[["x1", "x2", "x3", "x4"], [1.0, 1.0, 1.0, 1.0]],
            ]

    prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)


    prob.solve()


def test_CplexModel():
    '''using PyCPX'''
    from PyCPX import CPlexModel

    A = np.array([[1,0,0], [1,1,0], [1,1,1]])
    b = np.array([1,2,3])
    m = CPlexModel()

    x = m.new(3)
    t = m.new()

    m.constrain( abs((A*x - b)) <= t)
    m.minimize(t)
    m[x]

def test_linprog():
    import scipy
    from scipy.optimize import linprog
    c = [20, 30, 50] #Meksdmim
    A = [[-10, -20, -5], [-20, -15 ,-40]] #constrints
    b = [-55, -70]
    x0_bounds = (0, None)
    x1_bounds = (0, None)
    x2_bounds = (0, None)
    res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds, x2_bounds), options={"disp": True})
    return res


def Crispys_CplexModel(candidate_list, genes_list, df ,res_size):
    '''
    :param candidate_list:
    :param df: genes efficacy function
    :param res_size: number
    :param genes_list: list of lists. each sub list is a list of targets
    :return:
    '''
    #g_e_f = lambda targets_list : 1 - reduce(lambda x,y: x*y, map(lambda t: df(t)), targets_list) #genes efficacy function
    #make the M matrix
    #M = np.fromfunction(lambda i, j: 1 - reduce(lambda x,y: x*y, map(lambda t: df(candidate_list[i].seq, t), genes_list[j])), (len(candidate_list), len(genes_list)), dtype=float)
    #test np.fromfunction
    l1 = [5,6,7]
    l2 = [5,6,7]
    i,j= 0,0
    #map(lambda t: i+t, j)
    #M = np.fromfunction(lambda i, j: reduce(lambda x,y: x*y, map(lambda t: t, j)), (3,3), dtype=float)
    M = np.fromfunction(lambda i, j: list(map(lambda t: t+ i, j)), (3,3), dtype=float)
    #M = np.zeros((3,3))
    #M[(i,j)] = l1[i] + l2[j]
    print(M)

    ##or, fill M in a naive way:
    M = np.zeros(len(candidate_list), len(genes_list))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            fi = 1
            for target in genes_list[j]:
                fi = fi*(1-df(candidate_list[i].seq, target)) #df gives us the cut prob for each target
        M[i][j] = 1 - fi #the probability of the j gene to be cut by the i sgRNA. In the article, the sgRNA is j, and the gene is i.

    X = new(len(candidate_list), bool) #the candidates. if x[i] == 1, the i'th candidate is chosen
    m.constrain(X.sum()) <= res_size

    #Version 1:
    for j in range(len(t)): # for each gene
        t[j] = 1- reduce(lambda x,y: x*y,map(lambda i: 1 - M[i][j]*X[i], M.shape[0]))  #there is a problem here - it will be avaluated to a number, not stay with the variables. maybe a simpler way of using lamda experssion will work?

    m.maximise(t.sum())
    #version 2: maximise avg num of cuts, not of genes beaing cut:

    #t[i*j] = M[i][j]*X[i]

    #or, present this as a vectorial multipication:
    t = X.T() *  M
    m.maximise(t.sum())


def testCplexParliminaries():
    genes_list = [["AAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA"]] #one gene for now
    candidates_list = [Candidate("AAAAAAAAAAAAAAAAAAAA"), Candidate("AAAAAAAAAAAAAAAAAAAA")]
    Crispys_CplexModel(candidates_list, genes_list, cfd_funct, 1)


if __name__ == "__main__":
    #testCplexParliminaries()
    #r = test_linprog()
    #print(r)
    #test_CplexModel()
    test_Cplex()