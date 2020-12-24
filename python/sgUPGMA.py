__author__ = 'ItayM3'

##the old algirithm: buttoms up
import UPGMA
#hash table of the leaves

#the old algorithem


class sgGroup:
 def __init__(self, sgDict, genesList, seq):
  self.sgDict = sgDict # keys are genes names, values are lists of sequences
  self.genesList = genesList #might be more then one sgRNA per gene, so it's more readable to use the genesList as separated list
  self.seq = seq

 def __len__(self):
  return self.genesList.len

 def __str__(self):
  return "seq = " + self.seq+ "; genesList = " + str(self.genesList)

 def remove(self, elem):
  '''remove elem from the group'''
  elem_seq, elem_gene_name, elem_sgName = get_seq_and_gene_name_and_sg_name(elem)
  self.genesList.remove(elem_gene_name)
  del self.sgDict[elem_gene_name]


##getting the leaves
def test1():
 a = "aret"
 b = "ardw"
 c = "brdw"
 seq_list = [a,b,c]
 names = ["a", "b", "c"]
 '''
 matrix = UPGMA.make_initiale_matrix(UPGMA.d_f2,seq_list)
 #print(matrix)
 m2 = UPGMA.make_distance_matrix(names, matrix)
 #print(m2)
 m3 = m2.__repr__()
 #print(m3)
 upgma1 = UPGMA.make_UPGMA(m2)
 return upgma1
 '''
 return return_upgma(seq_list, names, UPGMA.d_f2)

def return_upgma(seq_list, names_list, df):
 '''input:  a list of names and a list of sequences, calibrated
  output: an upgma instance.
   '''
 matrix = UPGMA.make_initiale_matrix(df,seq_list)
 m2 = UPGMA.make_distance_matrix(names_list, matrix)
 #m3 = m2.__repr__()
 upgma1 = UPGMA.make_UPGMA(m2)
 return upgma1

def get_parent(root,clade):
 return root.get_path(clade)[-1]

def get_seq_and_gene_name_and_sg_name(elem):
 lst = elem.name.split(' ')
 gene_name = lst[0].split(':')[1]
 seq = lst[2].split(':')[1]
 sgName = lst[1].split(':')[1]
 return seq, gene_name, sgName

def make_singleton(elem):
 seq, gene_name, sgName = get_seq_and_gene_name_and_sg_name(elem)
 group = sgGroup({gene_name: seq},[gene_name],seq)
 return group


def add_singleton(elem, groups, names_hash_table):
 '''make a group with a single element'''
 group = make_singleton(elem)
 groups += [group]
 names_hash_table.remove(elem.name)

def add_to_group(elem, sgDict, genesList,names_hash_table):
 sgSeq, sgGeneName, sgName = get_seq_and_gene_name_and_sg_name(elem)
 if sgName in names_hash_table:
  if sgGeneName in sgDict:
   key = genesList[sgGeneName]
   key += [sgGeneName]
   sgDict[sgGeneName] = key
  else:
   sgDict[sgGeneName] = [sgSeq]
   genesList += [sgGeneName]
   names_hash_table.remove(sgName)

def DFS_like(elem, sgDict, genesList, names_hash_table, distance_from_leaf, delta, black_inner_nodes):
 distance_from_leaf += elem.branch_length
 if elem.name in black_inner_nodes: # this elem is not close enough to any other elem. nothing "good" will come out from it.
  return
 if distance_from_leaf >= delta/2: # this elem is not close enough to any other elem
  black_inner_nodes.add(elem.name)
  return
 if elem.is_terminal():
  add_to_group(elem, sgDict, genesList, names_hash_table)
 else:
  for a_clade in elem.clades:
   DFS_like(a_clade, sgDict, genesList, names_hash_table, distance_from_leaf, delta)



def theMain(delta, namesList, seqsList, df):
 '''the algorithem. df = distance function'''
 names_hash_table = set(namesList)# remove the leaf from here after finding him a group
 black_inner_nodes = set()
 groups = []
 #up1 = return_upgma(seq_list, names_list, df)
 up1 = test1() # a UPGMA tree
 iterator1 = up1.find_clades() #inisiating the chain
 elem = iterator1.__next__()
 while (not(elem.is_terminal())): #getting to a leaf
  elem = iterator1.__next__()
 if elem.branch_length >= (delta/2): # incloud '=' couse the distance function is a-symetric
  #group is a singletone
  add_singleton(elem, groups, names_hash_table)
  '''
  group = make_singleton(elem)
  groups += [group]
  names_hash_table.remove(elem.name)'''
 else:
  sgDict, genesList = {},[]
  distance_from_leaf= elem.branch_length
  while (distance_from_leaf < (delta/2)): # way delta/2? couse after getting to the parent, need to go back to the leaf, and it's UPGMA, so the distance to the leaves is symmetrical.
    if elem.is_terminal():
     add_to_group(elem, sgDict, genesList, names_hash_table)
    elem = get_parent(up1,elem)
    distance_from_leaf = distance_from_leaf + elem.branch_length
    DFS_like(elem, sgDict, genesList, names_hash_table, distance_from_leaf, delta/2, black_inner_nodes)
    black_inner_nodes.add(elem) #gone over all his children. IT'S TRUE THAT IN THIS POINT THIS NODE SHOULD TURN BLACK BECAUSE IT'S AN ULTRAMETRIC TREE
    add_singleton(elem,groups,names_hash_table)


#left to do: find the seq for each group, and return to the names_set the leaves that are not taken by the seq.
# how to find close sequences that are at opposite directions? is the Cas9 sensible to the direction- in that cas, nothing need to be done. how does global alingment reflect to it?
#done? find a more efficient way the go over the inner nodes, so an inner node which have no naked leaves left, wont be visited all the way to the leaves, but updating the "inner node colour to black" only after the leaves had been checked to be covered by the implementation. a way doing so: making a list of all the "black" inner nodes: which all their subtree is "black"

 #going over all the groups, choose a sequence for each, and for every sgRNA that isn't cough by this seq, make a singleton
 for group in groups:
  group_seq = choose_seq(group)
  group.seq = group_seq
  for elem in groups:
   if not(is_caught(group_seq,elem, df, delta)):
    make_singleton(elem)
    group.remove(elem)
 return groups

def choose_seq(group, naive = True):
 '''return the choosen seq. the one seq to capture them most'''
 if(naive):
  ''' do multiple alignment and find similar to consensus, in a navie way'''
     #make a suq list:
     sequences_for_alingment = []
     for gene in group.genesList:
      for current_seq in sgDict[gene]:
       sequences_for_alingment += [SeqRecord(Seq(current_seq, generic_dna), id= (gene + "seq: " +  current_seq))]
    align = MultipleSeqAlignment(sequences_for_alingment)
    #from here, can be done by loops, but, maybe there is a faster way, using saipy.

 return

def is_caught(group_seq, elem, df, delta):
 seq, gene_name, sgName = get_seq_and_gene_name_and_sg_name(elem)
 return df(group_seq,seq) < delta






#     while not a_clade.is_terminal(): #it's a binary tree with constant-clock assumption. no need for more testings
#      current_len_from__leaf = a_clade.branch_length
#      a_clade = a_clade.clades[i]
 #      add_to_group(elem, sgDict, genesList,names_hash_table) # that method also checks if this elem have already been removed from the names' set





#else: #need to get up at the tree, to the parent, and find all the node siblings
 # node_parent = get_parent(up1, elem)
  #now, cheack all the parent's children
if __name__ == "__main__":
 test11()