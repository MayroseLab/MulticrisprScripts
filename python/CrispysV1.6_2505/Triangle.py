import UPGMA
import random
import Metric

def is_pair_symetric(t1, t2, df = UPGMA.cfd_func):
	return  df(t1, t2) == df(t2, t1)

def test_symmery(df, N = 2000):
	Nuc = ["A", "G", "C", "T"]
	for i in range (N):
		t1 = ''.join(random.choice(Nuc) for _ in range(20))
		t2 = ''.join(random.choice(Nuc) for _ in range(20))
		if not is_pair_symetric(t1, t2, df):
			print("the function is not symmetric")
	print("the function is symmetric")


def test_triangle(df, N = 20000):
	Nuc = ["A", "G", "C", "T"]
	for i in range (N):
		t1 = ''.join(random.choice(Nuc) for _ in range(20))
		t2 = ''.join(random.choice(Nuc) for _ in range(20))
		t3 = ''.join(random.choice(Nuc) for _ in range(20))
		if not is_triangle_obteined(t1, t2, t3, df):
			print("the function is not triangle")
	print("the function obteined the triangle inequality")



def is_triangle_obteined(t1, t2, t3, df = UPGMA.cfd_func):
	d1, d2, d3 =  df(t1,t2), df(t2,t3), df(t1,t3)
	if (d1 + d2 < d3):
		print("d1: ", d1, "; d2: ", d2, "; d3: ", d3)
		return False
	return True

def test_triangle_by_pos(df, letter = 'A'):
	Nuc = ["A", "G", "C", "T"]
	t1, t2, t3 = [letter for i in range(20)], [letter for i in range(20)], [letter for i in range(20)]
	for i in range(20):
		for j in range(len(Nuc)):
			t1[i] = Nuc[(j+1)%4]
			t2[i] = Nuc[(j+2)%4]
			t3[i] = Nuc[(j+3)%4]
			#if not is_triangle_obteined(''.join(t1), ''.join(t2), ''.join(t3), df):
			if not is_triangle_obteined_after_S(''.join(t1), ''.join(t2), ''.join(t3), df):

				print("the function is not triangle")
				print(''.join(t1), ''.join(t2), ''.join(t3))
			t1[i], t2[i], t3[i] = letter, letter, letter
	print("if no print up to here, the function obteined the triangle inequality")



def test_symmetry_by_pairs(df, letter = 'A'):
	Nuc = ["A", "G", "C", "T"]
	t1, t2, t3 = [letter for i in range(20)], [letter for i in range(20)], [letter for i in range(20)]
	for i in range(20):
		for j in range(len(Nuc)):
			t1[i] = Nuc[(j+1)%4]
			t2[i] = Nuc[(j+2)%4]
			if not is_pair_symetric(''.join(t1), ''.join(t2), df):
				print("the function is not symmetric")
				print(df(''.join(t1), ''.join(t2)))
				print(df(''.join(t2), ''.join(t1)))
			t1[i] = Nuc[(j+3)%4]
			t2[i] = Nuc[(j)]
			if not is_pair_symetric(''.join(t1), ''.join(t2), df):
				print("the function is not symmetric")
				print(df(''.join(t1), ''.join(t2)))
				print(df(''.join(t2), ''.join(t1)))
	print("the function is symmetric")



def is_triangle_obteined_after_S(t1, t2, t3, df):
	d1_1, d1_2, d2_1, d2_2, d3_1, d3_2 =  df(t1,t2), df(t2,t1), df(t2,t3), df(t3, t2), df(t1,t3), df(t3,t1)
	d1, d2, d3 = (d1_1 + d1_2)/2,(d2_1 + d2_2)/2,(d3_1 + d3_2)/2
	if (d1 + d2 < d3):
		print("d1: ", d1, "; d2: ", d2, "; d3: ", d3)
		return False
	return True


if __name__ == "__main__":
	#test_symmery(UPGMA.shalem_score)
	#test_triangle(UPGMA.cfd_func)
	#test_triangle_by_pos(UPGMA.cfd_func)
	#test_symmetry_by_pairs(UPGMA.cfd_func)
    test_triangle_by_pos(Metric.find_dist_t)
