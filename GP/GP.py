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
		pi[node] = None
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
	i = 0
	while i < len(paths):
		while pi[paths[i][0]] != None:							# while 1st entry in path is not s
			path = paths.pop(i)									# remove that path
			paths.append([pi[path[0]][-1]] + path)				# Add back the path with the last predecessor in pi as new 1st entry
			for j in range(len(pi[path[0]])-1):					# For each predecessor other than the last one in the list:
				paths.append(deepcopy(path))					# Add a new copy of the path with that predecessor as new 1st entry
				paths[-1] = [pi[path[0]][j]] + paths[-1]
		i += 1
	if paths != [[t]]:
		return paths
	else:
		return None

# Implementation of Yen's K shortest paths based on Wikipedia pseudocode
# Inputs: G=graph, s=source node, t=target node, K=number of shortest paths
def yenKSP(Graph,s,t,K):
	G = deepcopy(Graph)
	k_paths = []
	potentials = []
	D,pi = dijkstra_all(G,s)
	paths = get_paths(pi,t)
	if paths is None:
		return None
	for path in paths:
		k_paths.append(path)
	start = 0
	stop = K-1
	for i in range(len(k_paths)):
		start += 1
		stop += 1
	for k in range(start,stop):						# To account for having 2+ tied paths from dijkstra_all
		for i in range(len(k_paths[k-1])-1):
			spur_node = k_paths[k-1][i]
			root_path = k_paths[k-1][:i+1]
			nodes_to_delete = []
			removed_edges = []
			for path in k_paths:
				if root_path == path[:i+1]:
					removed_edges.append((path[i],path[i+1],G[path[i]][path[i+1]]))
			for edge in removed_edges:
				if edge[1] in G[edge[0]].keys():
					del G[edge[0]][edge[1]]
				if edge[0] in G[edge[1]].keys():
					del G[edge[1]][edge[0]]
			for node in root_path:
				if node != spur_node:
					nodes_to_delete.append(node)
			deleted_nodes = del_node(G,nodes_to_delete)
			D_spur,pi_spur = dijkstra_all(G,spur_node)
			spur_paths = get_paths(pi_spur,t)
			total_path = deepcopy(root_path)
			if spur_paths is not None:
				for p in spur_paths:
					for node in p:
						if node not in total_path:
							total_path.append(node)
				if total_path not in potentials:
					potentials.append(total_path)
			for node,neighbors in deleted_nodes.items():
				G[node] = neighbors
				for v,w in neighbors.items():
					if v not in G.keys():
						G[v] = {}
					G[v][node] = w
			for edge in removed_edges:
				G[edge[0]][edge[1]] = edge[2]
				G[edge[1]][edge[0]] = edge[2]
		if len(potentials) > 0:
			potentials.sort()
			if potentials[0] not in k_paths:
				k_paths.append(potentials[0])
				potentials.pop()
	return k_paths

# Primary function - takes the graph, pos/neg labels for proteins, and a starting/ending node
# Returns a list of candidate fog pathway proteins
def get_candidates(Graph,L,s,t,K):
	G = deepcopy(Graph)
	candidates = []
	to_delete = []
	no_unknowns = False
	while not(no_unknowns) and len(candidates) <= 10:
		K_shortest_paths = yenKSP(G,s,t,K)
		if K_shortest_paths is None:
			break
		scores = []
		for p in K_shortest_paths:
			score = 0
			for node in p[1:len(p)-1]:			# Not counting fog/sqh
				if node in L.keys():
					if L[node] == "Positive":
						score += 1
			score = score / len(p[1:len(p)-1])
			scores.append(score)
		best_i = scores.index(max(scores))
		best_path = K_shortest_paths[best_i]
		for node in best_path[1:len(p)-1]:		# Not counting fog/sqh
			if node not in L.keys():
				candidates.append(node)
				to_delete.append(node)
		if len(to_delete) > 0:
			del_node(G,to_delete)
			to_delete = []
			print("Have "+str(len(candidates))+" candidate(s)")
		else:
			no_unknowns = True
	print("Final list of candidates: "+str(candidates))
	with open('Maddy_candidates.txt','w') as file:
	    for c in candidates[:len(candidates)-1]:
	    	file.write(c+"\n")
	    file.write(candidates[-1])

# Helper function for ease of use; makes a list of key-value pairs from the graph and then deletes
# those nodes and any edges containing those nodes from the graph
def del_node(G,node_list):
	deleted = {}
	for node in node_list:
		deleted[node] = deepcopy(G[node])
	for node in node_list:
		for v in G[node].keys():
			if v in G.keys():
				del G[v][node]
		del G[node]
	return deleted

# For some reason splits nodes into individual characters?????
# FIX THIS
def main():
	flyG,flyL = read_fly_interactome("interactome-flybase-collapsed-weighted.txt","labeled_nodes.txt")
	toyG,toyL = read_fly_interactome("toy_dataset.txt","toy_labeled.txt")
	toys = 'A1'
	toyt = 'G1'
	s = 'sqh'	# source node
	t = 'fog'	# target node
	K = 3		# number of shortest paths from s to t
	get_candidates(toyG,toyL,toys,toyt,K)
	#paths = yenKSP(flyG,s,t,K)
	#print("Final "+str(K)+" shortest paths:")
	#for p in paths:
	#	print(p)

# Results for fly interactome from desktop:
# 1. ['sqh', 'flw', 'Cul3', 'rdx', 'ci', 'sgg', 'Axn', 'dsh', 'Rho1', 'cta', 'fog']
# 2. ['sqh', 'Pp1-87B', 'flw', 'Cul3', 'rdx', 'ci', 'sgg', 'Axn', 'dsh', 'Rho1', 'cta', 'fog']
# 3. ['sqh', 'Pp1-87B', 'Mbs', 'flw', 'Cul3', 'rdx', 'ci', 'sgg', 'Axn', 'dsh', 'Rho1', 'cta', 'fog']

if __name__ == "__main__":
	main()