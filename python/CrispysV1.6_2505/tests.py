#import Covers
#import Covers1
import timeit
import os
import pickle
import sys


def run_exact_for_all():

	import shutil
	import sys, os
	sys.path.append("/groups/itay_mayrose/shiranabad/CRISPR/code/")
	from utilities import createJobFile


	#OUTPUT = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_5_11_PS/"
	path = "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/analysis_output_CDF_score_V1_5_11_PS"
	for dir in os.listdir(path):
		filepath = "/groups/itay_mayrose/shiranabad/CRISPR/gal/families_fasta/" + dir
		destdir = filepath #path + filename[:-3] + "/"
		dest_fasta = filepath#destdir + filename
		job_filename = createJobFile.create_job_file(dir, "python /groups/itay_mayrose/galhyams/CrispysV1.6/tests2.py "+ path + "/" + dir, dir + ".sh", "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/err_CDF_score_V1_5_11_PS/", "/groups/itay_mayrose/galhyams/CRISPIS_V_1.2_Data/jobs_CDF_score_V1_5_11_PS/")
		os.system('qsub -p -1 ' + job_filename)
		#break

if __name__ == "__main__":
	#print(*sys.argv)
	#print(len(*sys.argv[0]))
	#if len(*sys.argv) == 8: #"tests.py"
	run_exact_for_all()
