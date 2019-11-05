## Lab4 Starter Code

## this function returns two different graphs; one used in lecture and one used in lab.
def get_graph(example):
	'''
	Gets the edge list, adjacency list, and adjacency matrix of 
	example graphs. 
	Input: a string, either 'lecture' or 'lab'
	Output: node list, edge list (list of lists), adjacency list (dictionary), and ajacency matrix (list of lists)
	'''

	## DEFINE EDGE AND NODE LISTS
	if example == 'lecture':
		edge_list = [['A','B'],['A','C'],['B','C'],['B','E'],['B','D'],['C','E'],['D','E'],['E','F'],
			['C','H'],['C','G'],['G','H'],['C','H'],['I','J'],['K','L'],['L','M']]
		node_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M']
	elif example == 'lab':
		edge_list = [['B','D'],['A','D'],['A','C'],['A','E'],['B','C'],['C','E'],['E','F'],['F','J'],
			['J','L'],['I','J'],['H','I'],['A','H'],['G','H'],['G','K'],['H','K'],['K','M']]
		node_list = ['A','B','C','D','E','F','G','H','I','J','K','L','M']
	else:
		sys.exit('ERROR: example graph can be one of "lecture" or "lab". Exiting.')

	## CREATE ADJACENCY LIST
	adj_list = {n:[] for n in node_list}  ## another way to initialize dictonaries
	for e in edge_list:
		adj_list[e[0]].append(e[1])
		adj_list[e[1]].append(e[0])

	## CREATE ADJACENCY MATRIX (not super efficient, but it's OK 
	## because the graphs are small)
	adj_mat = []
	for i in range(len(node_list)):
		adj_mat.append([0]*len(node_list)) ## another way to add a row of 0's
	# adj_mat is now a matrix of zeros. 

	# Add ones where appropriate.
	for i in range(len(node_list)):  ## for every row...
		for j in range(i,len(node_list)): ## for the upper-right triangle...
			## we need to check whether the [i,j] edge OR the [j,i] edge
			## are in the edge_list.  
			edge = [node_list[i],node_list[j]]
			rev_edge = [node_list[j],node_list[i]]
			## if EITHER edge is in the list, update both the [i,j] and [j,i] entries.
			if edge in edge_list or rev_edge in edge_list:
				adj_mat[i][j] = 1
				adj_mat[j][i] = 1

	## print stats
	print('Graph "%s" has %d nodes and %d edges' % (example, len(node_list),len(edge_list)))
	print('Nodes: ',','.join(node_list))
	print('Edges: ',','.join(['%s-%s' % (e[0],e[1]) for e in edge_list]))
	print('\nAdjacency List (key --> value):')
	for key in adj_list:
		print('%s --> %s' % (key,','.join(adj_list[key])))
	print('\nAdjacency Matrix:')
	for row in adj_mat:
		print(' '.join([str(i) for i in row]))
		
	return node_list, edge_list, adj_list, adj_mat