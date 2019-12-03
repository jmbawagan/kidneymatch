'''
Creates donor-patient test data. This program creates one file:
	
	<name>.csv
	- contains the matrix / bipartite graph scores of the DP pairs
	
	Format:
	,P1,P2,P3,...
	D1,scoreD1P1,scoreD1P2,scoreD1P3,...
	D2,scoreD2P1,scoreD2P2,scoreD2P3,...
	D3,scoreD3P1,scoreD3P2,scoreD3P3,...
	...

	Given the following sample donor-patient test data,
				Patients
				P1   P2   P3
			D1   6    7    8
	Donors	D2   8    6    7
			D3   7    8    6

	This will be the corresponding file format:
	
		,patient1,patient2,patient3
		donor1,5,6,7
		donor2,7,5,6
		donor3,6,7,5

Command line arguments:

-n				number of donor-patient pairs
-f, --filename	name of the output file
-r, --random 	fills the matrix with random values (0-10)
-x				original donor pairs are incompatible (zero)
-v, --verbose	increase output verbosity

No options produces a 5x5 matrix filled with zeroes.

e.g.
python create_test_data.py -n 6 -r -x -v

'''

import argparse
import numpy as np

#>>Set arguments
parser = argparse.ArgumentParser(description="Creates donor-patient test \
								data. No options produces a 5x5 matrix \
								filled with zeroes.")

parser.add_argument("-n", type=int,	help="the number of donor-patient pairs")
parser.add_argument("-f", "--filename", help="name of the output file")
parser.add_argument("-r", "--random", action="store_true",
					help="fills the matrix with random values")
parser.add_argument("-x", action="store_true",
					help="all original DP pairs are incompatible")
parser.add_argument("-v","--verbose", action="store_true",
					help="increase output verbosity")

args = parser.parse_args()

#>>Set default values
n = 5
filename = "test.csv"

#>>Change n and filename if able
if args.n != None:
	n = args.n
	if n < 3:
		print("Invalid value for n, (n < 3).")
		raise SystemExit
if args.filename != None:
	filename = args.filename #+ ".csv"

#>>Create 2D array/matrix
if args.random:
	dp_matrix = np.random.randint(0,11,(n,n))
	if args.x:
		np.fill_diagonal(dp_matrix,0)
else:
	dp_matrix = np.zeros((n,n), dtype=int)


#>>Print output file
# format:
# 	,patient1,patient2,patient3
# 	donor1,5,6,7
# 	donor2,7,5,6
# 	donor3,6,7,5

fh = open(filename,"w")

for i in range(n):							#top row
	fh.write(",patient"+str(i+1))		
fh.write("\n")
for i in range(n):							#next rows
	fh.write("donor"+str(i+1))
	for j in range(n):
		fh.write(","+str(dp_matrix[i,j]))
	fh.write("\n")

fh.close()

#>>Verbose
if args.verbose:
	print("The bipartite graph/matrix:")
	print(dp_matrix)
	if args.random:
		if args.x:
			print("Random values set with all pairs incompatible.")
		else:
			print("Random values set.")

	
print("Created file",filename)
