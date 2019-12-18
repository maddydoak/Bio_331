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
from datetime import datetime
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

# From HW6, written by me, with modifications for this particular use case
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

# Implementation of Yen's K shortest paths based on Wikipedia pseudocode
# Inputs: G=graph (dict of dicts, see above), s=source node, t=target node, K=number of shortest paths
# Returns: list of lists of K shortest paths from s to t
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
	for k in range(start,stop):									# To account for having 2+ tied paths from dijkstra_all
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
			deleted_nodes = del_nodes(G,nodes_to_delete)
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

# Primary function,
# Inputs: Graph = graph (dict of dicts), L = dict of pos/neg labels, s = starting node, t = ending node, K = number of shortest paths to examine
# Returns: list of candidate fog pathway proteins
def get_candidates(Graph,L,s,t,K):
	G = deepcopy(Graph)
	candidates = []
	to_delete = []
	best_paths = []
	while len(candidates) < 100:
		K_shortest_paths = yenKSP(G,s,t,K)
		if K_shortest_paths is None:
			break
		print("K shortest paths: ",K_shortest_paths)
		def get_best_path(L,K_shortest_paths):
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
		paths_not_empty = True
		while len(to_delete) == 0 and paths_not_empty:
			best_path = get_best_path(L,K_shortest_paths)
			best_paths.append(best_path)
			for node in best_path[1:len(best_path)-1]:			# Not counting fog/sqh
				if node not in L.keys():
					candidates.append(node)
					to_delete.append(node)
			if len(to_delete) > 0:
				del_nodes(G,to_delete)
				if len(candidates) == 1:
					print(datetime.now().hour,":",datetime.now().minute,":",datetime.now().second," - ","Have ",str(len(candidates))," candidate")
				else:
					print(datetime.now().hour,":",datetime.now().minute,":",datetime.now().second," - ","Have ",str(len(candidates))," candidates")
			else:
				K_shortest_paths.remove(best_path)
				if len(K_shortest_paths) == 0:
					print(datetime.now().hour,":",datetime.now().minute,":",datetime.now().second," - ","Out of paths! Starting over...")
					paths_not_empty = False
		to_delete = []
	print("Final list of candidates: "+str(candidates))
	with open('Maddy_candidates_all.txt','w') as file:
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
	flyG,flyL = read_fly_interactome("interactome-flybase-collapsed-weighted.txt","labeled_nodes.txt")
	print("done. Getting best shortest paths...")
	s = 'sqh'	# source node
	t = 'fog'	# target node
	K = 5		# number of shortest paths from s to t
	best_paths = get_candidates(flyG,flyL,s,t,K)
	#graph_best_paths(graphspace,best_paths,flyL,flyG)

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
