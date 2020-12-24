## updated 19/10, 15:30
from __future__ import print_function
from Bio import SeqIO
import os
import heapq
import sys

####### for the new tool ###

def runIt(fullpath):
	callPerlOnGenes(fullpath)
	
def test_number_of_genes(fullpat,number):
	callPerlOnGenes(fullpath, number)
	
def testMultipleGenes(RNAinput,fullpath):
	##RNAFiles(RNAinput) # Done :)
##	fullpath= os.path.abspath(path)
	
	#dicti = part2(fullpath,2) #to be called later, on the other file

	writeDictiToFile(dicti)

def callPerlOnGenes(fullpath,number = sys.maxsize):
	i = 0
	for file in os.listdir(fullpath):
		if i>number:
			break
		if (file[-11:-1]== "RNAfile.tx"): ## a gene. hence: for each gene
			#if not os.path.exists(fullpath + '/' + file + "-holeGenomeWithExtras-s2n100pC"): #havan't been called yet
			i+=1
			callPerl(file)

def test_call(fullpath, number):
	i = 0
	for file in os.listdir(fullpath):
		
		if i>number:
			break
		if (file[-11:-1]== "RNAfile.tx"): ## a gene. hence: for each gene
			#if not os.path.exists(fullpath+'/' + file + "-holeGenomeWithExtras-s2n100pC"): #havan't been called yet
			i+=1
			print(file)


def callPerlOlder(gene):
	os.system("perl casot.pl -m=target -t="+gene+" -g=holeGenomeWithExtras.fa -o=tab -s=2 -n=100 -p=C -r=yes -l=19-20")

def callPerl(gene):
	#make files, including sh file
	print(gene)
	try:
		commandName = str(gene) + "command.txt"
		file = open(str(gene) + "command.txt",'w')
		file.write("perl /groups/itay_mayrose/galhyams/CasOT-1.0/casot.pl -m=target -t="+ str(gene) +" -g=holeGenomeWithExtras.fa -o=tab -s=2 -n=100 -p=N -r=no -l=20-20	callPerl"+ gene)
		file.close()
	except:
		print("failed at" + str(gene) + "command.txt" )
		sys.exit(-1)
	
	#call os.system
	string = "perl /groups/itay_mayrose/galhyams/CasOT-1.0/send_commands_to_queue_lecs_narrow.pl /groups/itay_mayrose/galhyams/CasOT-1.0/"+commandName+" /groups/itay_mayrose/galhyams/temp/ itaym"
	callQsub(string,gene) #don't know yet if the gene name is needed here
	
def callQsub(string,gene):
	os.system(string);	
	
	#system ("qsub -l itaim" +string+ ".sh");



def testcallPerl():
	os.system("perl casot.pl -m=target -t=Solyc00g005060.1.1RNAfile.txt -g=holeGenomeWithExtras.fa -o=tab -s=2 -n=100 -p=C -r=yes -l=19-20")



def writeDictiToFile(dicti):
	f = open("examlpe.txt",'w')
	for key,seq in dicti.items():
		f.write(str(key) + ' ' + str(seq) + '\n')
	f.close()

if __name__ == "__main__":
	fullpath = "/groups/itay_mayrose/galhyams/CasOT-1.0"
	test_number_of_genes(fullpath, 99999)
#	test_call(fullpath, 10)
	#callPerlOnAllGenes(fullpath)
	#test_number_of_genes(fullpath, sys.maxsize)
#	testMultipleGenes("4exampleGenes.txt","/groups/itay_mayrose/galhyams/CasOT-1.0")
