"""
Maddy Doak
"Group" Project
11/27/2019

Iteratively runs Dijkstra_all/get_paths from Fog to Sqh using costs (-log(weight)) of edges
Adds unlabeled nodes in the shortest path from Fog to Sqh to a list of candidates
Removes these unlabeled nodes from the graph at each iteration, then re-runs
Iteration ends when there are no more unlabeled nodes in the shortest path, OR the list has 15 candidates

OUTPUT: text file with list of candidates in order of rank
"""

from copy import deepcopy
from math import log

# From HW6, written by me
# Returns a dictionary of dictionaries of nodes and their neighbors/edge weights
# Also returns a dictionary of labeled nodes
def read_fly_interactome(graph_file,label_file):
	G = {}	# G = {u:{v:w,v2:w2,...},...}
	L = {}	# L = {u:Positive,v:Negative,...}
	with open(graph_file, 'rt') as file:
		next(file)	# Skip header line
		for line in file:
			line_list = line.strip().split()
			node1 = line_list[0]
			node2 = line_list[1]
			weight = float(line_list[2])
			cost = -log(weight,10)
			if node1 not in G:
				G[node1] = {node2:cost}
			else:
				G[node1][node2] = cost
			if node2 not in G:
				G[node2] = {node1:cost}
			else:
				G[node2][node1] = cost
	with open(label_file, 'rt') as f:
		next(f)
		for line in f:
			protein_info = line.strip().split()
			name = protein_info[0]
			label = protein_info[2]
			if name not in L.keys():
				L[name] = label
	return G,L

# Taken from HW6, written by myself
# Modified to take a graph as a dictionary of dictionaries
# Modified to fix bug - previously, distance was updated before checking whether the updated distance was less than OR equal to,
# so the case for the updated distance equaling the old distance was always used, never the case where updated < old
# Inputs: weighted graph G, a dictionary of dictionaries of nodes and their neighbors/edge weights {u:{v:w,v2:w2,...},...}
# and a starting node s
# Outputs: dictionaries D (distances between node pairs, {node:dist to node from s}) and pi (predecessors in shortest paths, {node:})
DEFAULT_DIST = 10000
def dijkstra_all(G,s):
	D = {}
	pi = {}
	for node in G.keys():
		D[node] = DEFAULT_DIST
		pi[node] = []
	D[s] = 0
	pi[s] = None
	to_visit = [s]
	while len(to_visit) != 0:
		min_node = None											# node with minimum distance value in D
		min_node_D = 100000
		for node in to_visit:
			if D[node] < min_node_D:
				min_node = node
				min_node_D = D[node]
		for neighbor in G[min_node].keys():						# For each neighbor of the node with the smallest distance from s
			if D[neighbor] == DEFAULT_DIST:						# If the neighbor hasn't been visited, add it to the "to_visit" list
				to_visit.append(neighbor)
			updated = float(D[min_node]) + G[min_node][neighbor]	# Possible updated distance from s to this neighbor of min_node = 
			if updated <= D[neighbor]:								# distance to min_node plus distance from min_node to this neighbor
				if updated < D[neighbor]:
					pi[neighbor] = [min_node]
				else:
					pi[neighbor].append(min_node)
				D[neighbor] = updated
		to_visit.remove(min_node)
	return D,pi

# From HW6, written by me
# pi = paths from s, t = target
# Returns a list of lists of nodes in the shortest path between the starting node s
# from dijkstra_all, and a target node t
def get_paths(pi,t):
	paths = [[t]]
	for i in range(len(paths)):
		while pi[paths[i][0]] != None:							# while 1st entry in path is not s
			path = paths.pop(i)									# remove that path
			paths.append([pi[path[0]][-1]] + path)				# Add back the path with the last predecessor in pi as new 1st entry
			for j in range(len(pi[path[0]])-1):					# For each predecessor other than the last one in the list:
				paths.append(deepcopy(path))					# Add a new copy of the path with that predecessor as new 1st entry
				paths[-1] = [pi[path[0]][j]] + paths[-1]
	return paths

# Implementation of Yen's K shortest paths based on Wikipedia pseudocode
# Inputs: G=graph, s=source node, t=target node, K=number of shortest paths
def yenKSP(Graph,s,t,K):
	G = deepcopy(Graph)
	k_paths = []
	potentials = []
	D,pi = dijkstra_all(G,s)
	paths = get_paths(pi,t)
	print(paths)
	for path in paths:
		k_paths.append(path)
	for k in range(5-len(k_paths)):						# To account for having 2+ tied paths from dijkstra_all
		for i in range(len(k_paths[k-1])-1):
			spur_node = k_paths[k-1][i]
			root_path = k_paths[k-1][:i+1]
			nodes_to_delete = []
			for path in k_paths:
				if root_path == path[:i+1]:
					nodes_to_delete.extend([path[i],path[i+1]])
			for node in root_path:
				if node != spur_node:
					nodes_to_delete.append(node)
			deleted_nodes = del_node(G,nodes_to_delete)
			D_spur,pi_spur = dijkstra_all(G,spur_node)
			spur_path = get_paths(pi_spur,t)
			total_path = root_path.extend(spur_path)
			potentials.append(total_path)
			for node,neighbors in deleted_nodes.items():
				G[node] = neighbors
				for v,w in neighbors.items():
					G[v][node] = w
		if len(potentials) != 0:
			potentials.sort()
			k_paths[k] = potentials[0]
			B.pop()
	return k_paths

# Helper function for ease of use
def del_node(G,node_list):
	deleted = {}
	for node in node_list:
		deleted[node] = G[node]
		for v in G[node].keys():
			del G[v][node]
		del G[node]
	return deleted

# Primary function - takes the graph, pos/neg labels for proteins, and a starting/ending node
# Returns a list of candidate fog pathway proteins
def get_candidates(Graph,L,s,t):
	G = deepcopy(Graph)
	candidates = []
	to_delete = []
	no_unknowns = False
	while not(no_unknowns) and len(candidates) <= 15:
		D,pi = dijkstra_all(G,s)
		paths = get_paths(pi,t)
		for path in paths:
			for node in path[1:len(path)-1]:
				if node not in L.keys():
					candidates.append(node)
					to_delete.append(node)
		if len(to_delete) > 0:
			for node in to_delete:
				for neighbor in G[node].keys():
					del G[neighbor][node]
				del G[node]
			to_delete = []
			print(candidates)
		else:
			no_unknowns = True
	with open('candidates.txt','w') as file:
	    for c in candidates[:len(candidates)-1]:
	    	file.write(c+"\n")
	    file.write(candidates[-1])

def main():
	G,L = read_fly_interactome("interactome-flybase-collapsed-weighted.txt","labeled_nodes.txt")
	s = 'sqh'	# source node
	t = 'fog'	# target node
	K = 3		# number of shortest paths from s to t
	#get_candidates(G,L,s,t)
	paths = yenKSP(G,s,t,3)

if __name__ == "__main__":
	main()