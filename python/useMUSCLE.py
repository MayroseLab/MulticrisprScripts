__author__ = 'ItayM5'
from Bio.Align.Applications import MuscleCommandline
import io
#from StringIO import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment

def call_MUSCLE1(fasta_input):
	'''input is a path to a file withe the sequnces at fasta format'''
	#help(MuscleCommandline)
	muscle_cline = MuscleCommandline(input = fasta_input)
	print(muscle_cline)
	stdout, stderr = muscle_cline()
	print(stdout)
	#align = AlignIO.read(io.StringIO(stdout), "fasta")
	#print(align)
	#return align

def call_MUSCLE2(input):
	'''input is a list of sequnces items'''

#records = (r for r in SeqIO.parse("opuntia.fasta", "fasta") if len(r) < 900)  ##insted of the input
	muscle_cline = MuscleCommandline(clwstrict=True)
	#muscle_cline = MuscleCommandline()
	#print(muscle_cline)
	#input = make_seqRecords(input)
	input = make_seqRecords2(input)
	handle = io.StringIO()
	SeqIO.write(input, handle, "fasta")
	data = handle.getvalue()

##You can then run the tool and parse the alignment as follows:

	stdout, stderr = muscle_cline(stdin=data)
	#align = AlignIO.read(io.StringIO(stdout))
	#print(align)
	#print(stdout)
	return stdout

def make_seqRecords(lst):
	'''input: lst is a list of tuples: (seq id, sequence)
	output: a list of seqRecord object, sutable to the input lst'''
	res = []
	for tup in lst:
		record = SeqRecord(Seq(tup[1],
                   IUPAC.IUPACAmbiguousDNA),
                   id=tup[0])
		res.append(record)
	return res

def make_seqRecords2(names_lst, seq_lst):
	'''input: 2 lst of the same length
	output: a list of seqRecord object, sutable to the input lst'''
	res = []
	for i in range(len(names_lst)):
		record = SeqRecord(Seq(seq_lst[i],
                   IUPAC.IUPACAmbiguousDNA),
                   id=names_lst[i])
		res.append(record)
	return res

def make_MSA_item(MSA_as_str):
	MSA_as_lst = MSA_as_str.split("\n")
	temp_MSA_lst = []
	for i in range(2, len(MSA_as_lst)-2):
		line = MSA_as_lst[i].split(" "*5)  ##number of " " : by the clwstrict format
		name = line[0]
		seq = line[1].strip()
		temp_MSA_lst.append(SeqRecord(Seq(seq, generic_dna), id=name))
	align = MultipleSeqAlignment(temp_MSA_lst,alphabet=generic_dna)
	print(align)
	return(align)


def return_MSA(input):
	'''
	:param input: a list of tuples. Each tuple is of the format (name, sequence)
	:return: A MSA item
	'''
	MSA_as_string = call_MUSCLE2(input)
	return make_MSA_item(MSA_as_string)


if __name__ == "__main__":
	sample1 = r"D:\Gal\MCdata\MUSCLE\test1.fa"
	sample2 = ["GTACGTTCAACGGGATCACTGGG","GCTTCACTAGAGAAGACTTCAGG","GGTGCTGATTGAAATGCTCCTGG","GTTGCGCTGAGGAAGAATTTGGG","GCCTGGAACGTGCCCACCAAGGG","GAACTTAGCACCAACTCTCCTGG","GGGACGATTACAGAATGAGAGGG"]
	sample3 = [("Solyc11g069040","GTACGTTCAACGGGATCACTGGG"), ("Solyc05g012430","GCTTCACTAGAGAAGACTTCAGG"),("Solyc06g008240","GGTGCTGATTGAAATGCTCCTGG"),("Solyc04g024690","GGGACGATTACAGAATGAGAGGG"),("Solyc06g066540","GTGGCTGTACCCACCGTGGTGGG")]
	sample4 = [("Solyc11g069040","GTACGTTCAACGGGATCACTGGG"),("Solyc05g012430","GCTTCACTAGAGAAGACTTCAGG"),("Solyc06g008240","GGTGCTGATTGAAATGCTCCTGG"),("Solyc01g110590","GTTGCGCTGAGGAAGAATTTGGG"),("Solyc12g037990","GCCTGGAACGTGCCCACCAAGGG")
				,("Solyc09g005530","GAACTTAGCACCAACTCTCCTGG"),("Solyc04g024690","GGGACGATTACAGAATGAGAGGG")]
	MSA = call_MUSCLE2(sample4)
	print(MSA)
	#call_MUSCLE1(sample1)
	make_MSA_item(MSA)