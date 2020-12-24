__author__ = 'GH'

import scipy
from scipy.stats.stats import pearsonr
import UPGMA

def fraseIt(inpath,df):
	'''
	:param inpath: the path of the exel file, in a txt tab dilimeter format
	:param df: distance fanction, or scoring function
	:return:
	'''
	sgVector = []
	sitesVector = []
	pamVector = []
	realScoreVector = []
	functionScoreVector = []
	f = open(inpath,"r")
	next(f)
	for line in f:
		line_as_array = line.split("\t")
		sgVector.append(line_as_array[2])
		sitesVector.append((line_as_array[3]))
		pamVector.append(line_as_array[4])  #might want to consider only NGG pam when calculating the pearson corilation ->removing unsutable sequnces from the vectors.
		realScoreVector.append(float(line_as_array[11]))
		functionScore = df(line_as_array[2], line_as_array[3])
		functionScoreVector.append(functionScore)
	##regresion time!!
	#pearson = scipy.stats.pearsonr(realScoreVector, functionScoreVector)
	pearson = pearsonr(realScoreVector, functionScoreVector)  ## returns a tuple. first argument is the pearson corelation. the second is 2-tailed p-value
	return pearson

if __name__ == '__main__':
	df1 = UPGMA.p_distance
	df2 = UPGMA.shalem_score
	folder = r"D:\Gal\MC evaluate scoring funcitions"
	inpath = folder + "\parsed_htgts_ots.txt"
	print(fraseIt(inpath, df1))
	#print(fraseIt(inpath, df2)) #only works when the str are alweys 20bp of length