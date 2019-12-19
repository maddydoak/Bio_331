"""
Maddy Doak
"Group" Project

DESCRIPTION:
Iteratively runs Yen's KSP from Fog to Sqh using costs (-log(weight)) of edges
Of the K shortest paths, picks the one with the highest proportion of positively-labeled nodes
Adds unlabeled nodes in this path to a list of candidates
Removes these unlabeled nodes from the graph, then re-runs
Iteration ends when there are no more unlabeled nodes in any of the K shortest paths, OR the list has >10 candidates

INSTRUCTIONS:
To start, go to the main() function at the bottom.
In the "read_fly_interactome" function, add filepaths/filenames for the .txt files with graph nodes/edges and node labels, respectively,
if different from those already entered (these may be found in the "inputs" Google Drive folder).
You may change the source/target nodes (s and t, respectively) if desired, but this script is meant to use sqh/fog as the source and target.
You may adjust K to be any number of paths you want to choose from at each iteration, but 5 is a good starting point.
Nothing else need be changed in order to run the program as you normally would.

OUTPUT: text file with list of candidates in order of rank
"""

from copy import deepcopy
from math import log
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

# From HW6, written by me, with modifications for this particular use case
# Returns a dictionary of dictionaries of nodes and their neighbors/edge weights
# Also returns a dictionary of labeled nodes
def read_fly_interactome(graph_file,label_file):
	G = {}	# G = {u:{v:w,v2:w2,...},...}
	L = {}	# L = {u:Positive,v:Negative,...}
	with open(graph_file, 'rt') as file:						# Get node/edge info for graph
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
	with open(label_file, 'rt') as f:							# Get node labels
		next(f)
		for line in f:
			protein_info = line.strip().split()
			name = protein_info[0]
			label = protein_info[2]
			if name not in L.keys():
				L[name] = label
	return G,L

# Taken from HW6, written by myself; modified to take a graph as a dict of dicts
# Modified to fix bug - before, distance was updated before checking if updated distance was < OR =, so always used = case
# Inputs: G (weighted/undirected graph as a dict of dicts of nodes + their neighbors/edge weights {u:{v:w,v2:w2,...},...}), and starting node s
# Outputs: dicts D (distances between node pairs, {node:dist from s}) and pi (predecessors in shortest paths, {node:[pred,alt_pred]})
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
			updated = float(D[min_node]) + G[min_node][neighbor]# Possible updated distance from s to this neighbor of min_node =
			if updated <= D[neighbor]:							# distance to min_node plus distance from min_node to this neighbor
				if updated < D[neighbor]:
					pi[neighbor] = [min_node]
				else:
					pi[neighbor].append(min_node)
				D[neighbor] = updated
		to_visit.remove(min_node)
	return D,pi

# From HW6, written by me
# Inputs: pi = paths from s, t = target node
# Returns: list of lists of nodes in shortest path from starting node s (from dijkstra_all) and target node t
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

# Implementation of Yen's K shortest paths based on Wikipedia pseudocode; see Wiki page for description of the "spur node," "root path," etc
# Inputs: G=graph (dict of dicts, see above), s=source node, t=target node, K=number of shortest paths
# Returns: list of lists of K shortest paths from s to t
def yenKSP(Graph,s,t,K):
	G = deepcopy(Graph)
	k_paths = []												# List to be returned with the K shortest paths from s to t
	potentials = []												# Possible shortest paths to be added to the above list
	D,pi = dijkstra_all(G,s)
	paths = get_paths(pi,t)										# Function that takes distances and predecessors from Dijkstra and returns a list of lists of paths
	if paths is None:
		return None												# For the case of an unconnected graph; exits yenKSP
	for path in paths:
		k_paths.append(path)									# Starts by adding the shortest path (or paths if tied) from Dijstra to the main list
	start = 0
	stop = K-1
	for i in range(len(k_paths)):								# Adjusts the range of the main for loop to correct for possible ties from Dijkstra
		start += 1
		stop += 1
	for k in range(start,stop):									# The main iterative loop that looks for the next-shortest paths
		for i in range(len(k_paths[k-1])-1):					# Looks at the last shortest path in the main list, except for the last node
			spur_node = k_paths[k-1][i]							# a node in the current path to try and "branch off from" and create a new next-shortest path
			root_path = k_paths[k-1][:i+1]						# the path up to and including the spur node
			nodes_to_delete = []								# Passed to function that deletes the node and all edges containing it from the graph
			removed_edges = []									# List of edges that were temporarily removed to be added back later
			for path in k_paths:
				if root_path == path[:i+1]:						# Checks all other paths in main list to see if they have the same "root"
					removed_edges.append((path[i],path[i+1],G[path[i]][path[i+1]]))		# If so, removes the edge connecting the root to the rest of the path
			for edge in removed_edges:
				if edge[1] in G[edge[0]].keys():
					del G[edge[0]][edge[1]]
				if edge[0] in G[edge[1]].keys():
					del G[edge[1]][edge[0]]
			for node in root_path:								# Removes all non-spur nodes in the root path from the graph
				if node != spur_node:
					nodes_to_delete.append(node)
			deleted_nodes = del_nodes(G,nodes_to_delete)
			D_spur,pi_spur = dijkstra_all(G,spur_node)			# Runs Dijkstra using the spur node as the source
			spur_paths = get_paths(pi_spur,t)					# And gets the shortest path from the spur to the target
			total_path = deepcopy(root_path)
			if spur_paths is not None:							# Creates a new potential shortest path with the root + the new branched off path from spur to t
				for p in spur_paths:
					for node in p:
						if node not in total_path:
							total_path.append(node)
				if total_path not in potentials:
					potentials.append(total_path)
			for node,neighbors in deleted_nodes.items():		# This and the next for loop add back the removed nodes/edges
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
			if potentials[0] not in k_paths:					# Sorts potentials by length, adds the shortest to k_paths
				k_paths.append(potentials[0])
				potentials.pop()
	return k_paths

# Primary function, runs Yen's KSP repeatedly and removes unlabeled nodes, adding them to a list of candidates
# Inputs: Graph = graph (dict of dicts), L = dict of pos/neg labels, s = starting node, t = ending node, K = number of shortest paths to examine
# Returns: list of all generated best paths (list of candidates is written to a .txt doc, not returned)
def get_candidates(Graph,L,s,t,K):
	G = deepcopy(Graph)
	candidates = []												# The list of unlabeled nodes removed from the best paths
	best_paths = []												# A list of the chosen best paths, returned so that they can be viewed if desired
	to_delete = []												# Keeps track of nodes to be removed from the graph at each iteration
	paths_not_empty = True										# Used in while loops, becomes False if there are no unlabeled nodes in any of the K shortest paths
	while len(candidates) < 10 and paths_not_empty:				# This value may be modified for a longer list of candidates; this is for a short list
		K_shortest_paths = yenKSP(G,s,t,K)
		if K_shortest_paths is None:
			break
		def get_best_path(L,K_shortest_paths):					# Picks the path out of the K shortest paths that has the highest proportion of positive-labeled nodes
			scores = []
			for p in K_shortest_paths:
				score = 0
				for node in p[1:len(p)-1]:						# Not counting fog/sqh
					if node in L.keys():
						if L[node] == "Positive":
							score += 1
				score = score / len(p[1:len(p)-1])
				scores.append(score)
			return K_shortest_paths[scores.index(max(scores))]
		while len(to_delete) == 0 and paths_not_empty:			# While still trying to find a path with unlabeled nodes...
			best_path = get_best_path(L,K_shortest_paths)
			best_paths.append(best_path)						# All best paths are tracked, regardless of if they generate candidates
			for node in best_path[1:len(best_path)-1]:			# Not counting fog/sqh
				if node not in L.keys():
					candidates.append(node)
					to_delete.append(node)						# If the node isn't in the L dictionary, that means it is unlabeled, and thus a candidate
			if len(to_delete) > 0:
				del_nodes(G,to_delete)
			else:												# If there are no unlabeled nodes, remove the best path from the list of K shortest paths
				K_shortest_paths.remove(best_path)
				if len(K_shortest_paths) == 0:					# If there are no paths left in the list of K shortest, ends this while loop and the main one
					paths_not_empty = False						# This could potentially be changed so that K is increased and the main loop is run again
		to_delete = []
	print("Final list of candidates: "+str(candidates))
	with open('Maddy_candidates_all.txt','w') as file:			# Writes the results to two files, one with the top 10, and one with all if there are >=10
		for c in candidates[:len(candidates)-1]:
			file.write(c+"\n")
		file.write(candidates[-1])
	if len(candidates) >=10:
		with open('Maddy_candidates_10.txt','w') as file:
			for c in candidates[:9]:
				file.write(c+"\n")
			file.write(candidates[9])
	return best_paths

# Helper function for ease of use; makes a list of key-value pairs from the graph and then deletes
# those nodes and any edges containing those nodes from the graph
def del_nodes(G,node_list):
	deleted = {}
	for node in node_list:
		deleted[node] = deepcopy(G[node])
	for node in node_list:
		for v in G[node].keys():
			if v in G.keys():
				del G[v][node]
		del G[node]
	return deleted

# Posts or updates existing graph on GraphSpace
def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

# Graphs network of shortest paths between s, t
def graph_best_paths(gs_session,paths,labels,graph):
	G = GSGraph()
	G.set_name("Doak - Best Paths from Sqh to Fog")
	G.set_tags(['GP'])
	st_nodes = [paths[0][0],paths[0][-1]]
	nodes = {}
	edges = []
	n_colors = ["#e76af7","#b767f5","#7666f2","#67aef5","#65d7eb","#62e391","#ade065"]
	for p in paths:
		for i in range(len(p)-1):
			if p[i] not in nodes.keys():
				nodes[p[i]] = paths.index(p)
			if (p[i],p[i+1]) not in edges and (p[i+1],p[i]) not in edges:
				edges.append((p[i],p[i+1],graph[p[i]][p[i+1]]))
		if p[-1] not in nodes.keys():
			nodes[p[-1]] = paths.index(p)
	for node,i in nodes.items():
		G.add_node(node,label=node)
		if node in st_nodes:
			G.add_node_style(node,
							 color = "#fa6bac",
							 height = 30,
							 width = max(30,15*(len(node))))
		else:
			if node in labels.keys():
				if labels[node] == "Positive":
					G.add_node_style(node,
									 color = n_colors[nodes[node]],
									 border_width = 3,
									 border_color = 'green',
									 style = 'dashed',
									 height = 30,
									 width = max(30,15*(len(node))))
				elif labels[node] == "Negative":
					G.add_node_style(node,
									 color = n_colors[nodes[node]],
									 border_width = 3,
									 border_color = 'red',
									 style = 'dashed',
									 height = 30,
									 width = max(30,15*(len(node))))
			else:
				G.add_node_style(node,
								 color = n_colors[nodes[node]],
								 border_width = 5,
								 height = 30,
								 width = max(30,15*(len(node))))
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],width=1+float(edge[2]))
	G.set_data(data={'description': 'shortest paths in fly interactome from sqh to fog, with sqh and fog highlighted in'\
	 					' pink. First shortest paths are warmer (closer to pink), shifting to purple, blue, then green'\
						' in order of when the paths were generated after removing previous nodes. Red-bordered nodes'\
						' are negative, green-bordered nodes are positive, and black-bordered nodes are unlabeled candidates.'})
	post(G, gs_session)

def main():
	print("Reading fly interactome...")
	G,L = read_fly_interactome("interactome-flybase-collapsed-weighted.txt","labeled_nodes.txt")
	print("done. Getting best shortest paths...")
	s = 'sqh'	# source node
	t = 'fog'	# target node
	K = 5		# number of shortest paths from s to t
	best_paths = get_candidates(G,L,s,t,K)

"""
K=3

Final candidates:
['flw', 'rdx', 'ci', 'Axn', 'Pp1-87B', 'l(2)gl', 'numb', 'N', 'CycK', 'Cdk2', 'Dsor1']

Final best paths:
['sqh', 'flw', 'Cul3', 'rdx', 'ci', 'sgg', 'Axn', 'dsh', 'Rho1', 'cta', 'fog']
['sqh', 'Pp1-87B', 'Mbs', 'Rho1', 'cta', 'fog']
['sqh', 'Rho1', 'cta', 'fog']
['sqh', 'zip', 'l(2)gl', 'par-6', 'aPKC', 'numb', 'N', 'dsh', 'Rho1', 'cta', 'fog']
['sqh', 'Rho1', 'cta', 'fog']
['sqh', 'zip', 'sgg', 'arm', 'dsh', 'Rho1', 'cta', 'fog']
['sqh', 'Rho1', 'cta', 'CycK', 'Cdk2', 'rl', 'Dsor1', 'Raf', 'tor', 'fog']
"""

"""
K=5

Final candidates:
['Pp1-87B', 'flw', 'Roc1b', 'ci', 'CkIalpha', 'Axn', 'l(2)gl', 'numb', 'N', 'CycK', 'Cdk2', 'Dsor1']

Final best paths:
[['sqh', 'Pp1-87B', 'Mbs', 'Rho1', 'cta', 'fog'],
['sqh', 'flw', 'Cul3', 'Roc1b', 'ci', 'CkIalpha', 'arm', 'Axn', 'dsh', 'Rho1', 'cta', 'fog'],
['sqh', 'Rho1', 'cta', 'fog'],
['sqh', 'zip', 'l(2)gl', 'par-6', 'aPKC', 'numb', 'N', 'dsh', 'Rho1', 'cta', 'fog'],
['sqh', 'Rho1', 'cta', 'fog'],
['sqh', 'zip', 'sgg', 'arm', 'dsh', 'Rho1', 'cta', 'fog'],
['sqh', 'Rho1', 'cta', 'CycK', 'Cdk2', 'rl', 'Dsor1', 'Raf', 'tor', 'fog']]

"""

if __name__ == "__main__":
	main()
