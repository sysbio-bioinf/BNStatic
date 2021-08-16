print("FVSsage.sage file running!")

#The dataset of 35 networks was previously iteratively reduced by removing sink and source nodes. These reduced networks are then used to calculate the FVS.

#import libraries
import csv
from numpy import array


o = open('/home/weidner/Desktop/SageScript/outputFVS-selfloopinputsremoved-1indexed.txt','w')

for i in range(35):
	print("i = " + str(i))
	n = i + 1 #i counts from 0 to 34, need 1 to 35 for nets
	#load txt/csv file nr n
	print("Network number " + str(n)) 
	file=r'/home/weidner/Desktop/SageScript/Reduced adjmats noSelfLoop/net_' + str(n) + '.csv'
	reader=csv.reader(open(file))
	L=[]
	for row in reader:                   
   	 L.append(row)

	GeneNames = str((L[0])[1:]) #first row of matrix should contain gene names
	print(type(GeneNames))

	#convert to Matrix
	M = matrix(RDF,L[1:])

	#print(type(M)) #<class 'sage.matrix.matrix_real_double_dense.Matrix_real_double_dense'>
		# Is integer matrix necessary for adjmat->graph conversion?
	#convert to integer matrix (sage.matrix.matrix_integer_dense.Matrix_integer_dense) ?

	M = M.delete_columns([0]) #slice off column 0 from matrix so that only adjmat remains
	#show(M)

	#convert reduced adjmat to DiGraph
	DiG = DiGraph(M, vertex_labels=GeneNames, loops = True)
	
	reducedSize = DiG.num_verts()
	FVS = sage.graphs.generic_graph.GenericGraph.feedback_vertex_set(DiG)
	FVS_1indexed = [ i+1 for i in FVS ]
	print("FVS, 1-indexed, indices of reducedAdjmat: " + str(FVS_1indexed))

	o.write(str("Network number " + str(n)))
	o.write("\n")
	o.write(str(DiG))
	o.write("\n")
	o.write(str(FVS_1indexed))
	o.write("\n\n")





