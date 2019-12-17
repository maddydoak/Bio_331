"""
Implementation of Yen's KSP
My graph structure is {node:{neighbor:edge_weight,neighbor2:edge_weight2,...},...}
where all nodes are included as keys in the graph, and all nodes have all of their neighbors in their dictionary of values
(So there are duplicates of edges)
"""

# Implementation of Yen's K shortest paths based on Wikipedia pseudocode
# Inputs: G=graph (dict of dicts), s=source node, t=target node, K=number of shortest paths
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
