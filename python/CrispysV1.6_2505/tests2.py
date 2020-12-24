import Covers2
#import Covers
import sys
if __name__ == "__main__":
	print(*sys.argv[1:])
	Covers2.run_exact_for_famliy(str(*sys.argv[1:]))
	#Covers.run_exact_for_famliy(str(*sys.argv[1:]))