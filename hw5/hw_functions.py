## Run Dijkstra's in the weighted, undirected graph.
## INPUT: set of nodes, 3-element list of edges [node1,node2,weight], source s
## OUTPUT: Dictionary of distances (D), Dictionary of predecessors (pi)
def dijkstra(nodes,edges,s):
	## Build adjacency list that contains the weights of the edge.
	## e.g., for edge (u,v), you can access the weight of that edge
	## with adj_list[u][v] OR adj_list[v][u]
	adj_list = get_adj_list_with_weights(edges)

	LARGE_NUM = 1000000 ## like "infinity" here.

	## initialize distances dictionary D.
	D = {n:LARGE_NUM for n in nodes}

	## initialize predecessor dictionary pi.
	pi = {n:None for n in nodes}

	## set distance to s to be 0
	D[s] = 0

	## Queue is a dictionary (slow implementation)
	## This could be sped up with a proper priority queue,
	## but is fine for this homework.
	## The queue values start as the distances for each node.
	Q = {n:D[n] for n in nodes} 
	
	while len(Q) > 0: ## While we haven't visited all the nodes...
		## Find the node with the minimum weight.
		w = None 
		for n in Q: ## for every node in the Queue...
			if w == None or Q[n] < Q[w]: ## if we haven't set w yet or n is better...
				w = n ## set w to be this node.

		## remove w from queue
		del Q[w] 
		
		## Iterate through the neighbors of w
		for x in adj_list[w]:
			## If the current distance to x is larger than coming from w, update
			if D[x] > D[w] + adj_list[w][x]:
				D[x] = D[w] + adj_list[w][x] ## update the distance
				pi[x] = w ## update the predecessor (we came from w)
				Q[x] = D[x] ## update the entry in the queue
				
	return D,pi

## Make an adjacency list that contains the weights of each edge.
## e.g., for edge (u,v), you can access the weight of that edge
## with adj_list[u][v] OR adj_list[v][u]
## INPUT: 3-element list of edges [node1,node2,weight]
## OUTPUT: dictionary of dictionaries
def get_adj_list_with_weights(edges):
	adj_list = {}
	## loop over all edges.
	for u,v,w in edges: ## another way to specify elements of key

		## We want to add the key-value pair (v,w) to adj_list[u].
		## First see if u is a key in adj_list.
		if u not in adj_list:
			adj_list[u] = {}  ## add the key (value is a DICTIONARY)
		## Add the key-value pair (v,w) to adj_list[u]
		adj_list[u][v] = w

		## We want to add the key-value pair (u,w) to adj_list[v].
		## First see if v is a key in adj_list.
		if v not in adj_list:
			adj_list[v] = {}  ## add the key (value is a DICTIONARY)
		## Add the key-value pair (u,w) to adj_list[v]
		adj_list[v][u] = w
		
	return adj_list


## Given a predecessor dictionary (e.g, pi from dijkstra()) and 
## a node, return the path as a list of nodes.
## If you copy this function, you MUST add comments that describes
## what it's doing.
def get_path(pi,node):
 	path = [node]
 	while pi[path[0]] != None:
 		path = [pi[path[0]]] + path
	return path
