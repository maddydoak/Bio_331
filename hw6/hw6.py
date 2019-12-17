"""
Maddy Doak
BIO 331 - HW6
Girvan-Newman algorithm for community detection on undirected, weighted graphs
"""
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
import copy

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

DEFAULT_DIST = 10000

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def post_gn_graph(graph,communities,k,name,gs_session,actual_groups = None):
	G = GSGraph()
	G.set_name(name)
	G.set_tags(['HW6'])
	for node in graph.keys():
		G.add_node(node,label=node)
		colors = ["#FE4365","#FC9D9A","#F9CDAD","#C8C8A9","#83AF9B","#82a4ad","#8c91b8","#997ca6"]
		shapes = ['ellipse','rectangle','triangle','diamond','star','hexagon','heptagon','octagon']
		comm = communities[k]
		node_id = -1
		for c in comm:
			if node in c:
				node_id = comm.index(c)
				break
		G.add_node_style(node,
						 color = colors[node_id],
						 height = 30,
						 width = 30)
		if actual_groups is not None:
			G.add_node_style(node,
							 shape = shapes[int(actual_groups[node])-1],
						 	 color = colors[node_id],
						 	 height = 30,
						 	 width = 30)
	edges = []
	for node,neighbors in graph.items():
		for neighbor in neighbors:
			if (node,neighbor[0],neighbor[1]) not in edges and (neighbor[0],node,neighbor[1]) not in edges:
				edges.append((node,neighbor[0],neighbor[1]))
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],width=1+float(edge[2]))
	G.set_data(data={'description': str(k)+' partitions using the Girvan-Newman algorithm on a weighted, undirected graph. Shapes for'\
									' actual badger graphs are: 1=ellipse, 2=rectangle, 3=triangle, 4=diamond, 5=star, 6=hexagon,'\
									' 7=heptagon, 8=octagon'})
	post(G, gs_session)

# Reads in the actual badger setts so that nodes can be assigned shapes
def read_badger_groups(filename):
	G = {}									# G = {u:[(v,w),(v2,w2)]}
	with open(filename, 'rt') as file:
		for line in file:
			badger = line.strip().split()
			badger_ID = badger[0]
			badger_group = badger[3]
			G[badger_ID] = badger_group
	return G

# Returns a dictionary (see comment below) of nodes and their neighbors/edge weights
def get_graph(filename):
	G = {}									# G = {u:[(v,w),(v2,w2)]}
	with open(filename, 'rt') as file:
		for line in file:
			edge = line.strip().split()
			node1 = edge[0]
			node2 = edge[1]
			weight = edge[2]
			if node1 not in G:
				G[node1] = [(node2,weight)]
			else:
				G[node1].append((node2,weight))
			if node2 not in G:
				G[node2] = [(node1,weight)]
			else:
				G[node2].append((node1,weight))
	return G

# Returns a dictionary of communities (see format below) based on the Girvan-Newman
# clustering algorithm, which removes the edges with highest edge betweenness and then
# determines whether new graph partitions have been created
# Prints each community as it is determined
def girvan_newman(Graph):
	G = copy.deepcopy(Graph)
	communities = {}						# {k:[partitions]}
	all_nodes = []
	for node in G.keys():
		all_nodes.append(node)
	communities[1] = [all_nodes]
	print(communities[1])
	previous = 1
	while len(communities) < len(all_nodes):
		eb_dict = edge_betweenness(G)		# {(u,v):EB (int)}
		edge_to_remove = [(-1,-1),-1]		# edge with the highest edge betweenness score, which is removed
		for edge,eb in eb_dict.items():
			if eb > edge_to_remove[1]:
				edge_to_remove = [edge,eb]
		node1 = edge_to_remove[0][0]
		node2 = edge_to_remove[0][1]
		if node1 != -1 and node2 != -1:
			for neighbor in G[node1]:
				if neighbor[0] == node2:
					G[node1].remove(neighbor)
			for neighbor in G[node2]:
				if neighbor[0] == node1:
					G[node2].remove(neighbor)
		D,pi = dijkstra_all(G,node1)
		if D[node2] == DEFAULT_DIST:		# If removing an edge split a partition (node1 can no longer reach node2)
			new1 = []						# Then two new partitions are made based on which nodes node1 can still reach (new1)
			new2 = []
			curr_comm = 0
			for community in communities[previous]:
				if node1 in community and node2 in community:
					curr_comm = community
					break
			for node,dist in D.items():
				if node in curr_comm:
					if dist != DEFAULT_DIST:
						new1.append(node)
					else:
						new2.append(node)
			communities[previous+1] = copy.deepcopy(communities[previous])
			communities[previous+1].remove(curr_comm)
			communities[previous+1].append(new1)
			communities[previous+1].append(new2)
			print(communities[previous+1])
			previous += 1
	return communities

# Calculates the edge betweenness score for each edge in the graph
def edge_betweenness(G):
	eb = {}									# {(u,v):eb_score (int)}
	ignoreu = []							# So duplicates (e.g. (u,v) and (v,u)) are not checked
	for u,neighbors in G.items():
		ignoreu.append(u)
		for vw in neighbors:
			v = vw[0]
			if v not in ignoreu:
				eb[(u,v)] = 0				# The actual EB score, with a value added for each path from s to t
				ignore = []					# Also to prevent duplicate paths being checked
				for s in G.keys():
					ignore.append(s)
					D,pi = dijkstra_all(G,s)
					for t in G.keys():
						if t not in ignore:
							st_paths = get_paths(pi,t)
							paths_through_uv = 0
							for path in st_paths:
								if u in path:
									if path.index(u) > 0:
										if path[path.index(u)-1] == v:
											paths_through_uv += 1
									if path.index(u) < len(path)-1:
										if path[path.index(u)+1] == v:
											paths_through_uv += 1
							eb[(u,v)] += paths_through_uv/len(st_paths)
	return eb

# Returns all possible shortest paths from a starting node s to all other nodes, including ties
def dijkstra_all(G,s):
	D = {}
	pi = {}
	for node in G.keys():
		D[node] = DEFAULT_DIST
		pi[node] = []
	D[s] = 0
	visited = [s]
	while len(visited) != 0:
		min_node = ['NA',100000]		# [node name, minimum distance value]
		for node in visited:
			if D[node] < min_node[1]:
				min_node = [node,D[node]]
		for neighbor in G[min_node[0]]:
			if D[neighbor[0]] == DEFAULT_DIST:
				visited.append(neighbor[0])
			updated = float(D[min_node[0]]) + float(neighbor[1])
			if D[neighbor[0]] >= updated:
				D[neighbor[0]] = updated
				if D[neighbor[0]] > updated:
					pi[neighbor[0]] = min_node[0]
				else:
					pi[neighbor[0]].append(min_node[0])
		visited.remove(min_node[0])
	return D,pi

# pi = paths from s, t = target
# Returns a list of lists of nodes in the shortest path between the starting node s
# from dijkstra_all, and a target node t
def get_paths(pi,t):
	paths = [[t]]
	for i in range(len(paths)):
		while len(pi[paths[i][0]]) > 0:
			path = paths.pop(i)
			paths.append([pi[path[0]][-1]] + path)
			for j in range(len(pi[path[0]])-1):
				paths.append(copy.deepcopy(path))
				paths[-1] = [pi[path[0]][j]] + paths[-1]
	return paths

def main():
	# Calculates the k = 4 partitions when partitioning the toy graph and posts the results
	# with each partition colored differently
	toy_graph = get_graph('toy_dataset.txt')
	communities = girvan_newman(toy_graph)
	post_gn_graph(toy_graph,communities,4,'Doak - HW6 Toy Graph',graphspace)

	# Calculates the k = 8 partitions for the badger graph, posting the graph with each
	# calculated partition being a different color, and the actual badger setts represented
	# by the node shapes
	badger_graph = get_graph('badger_edge_costs.txt')
	actual_badg = read_badger_groups('badger-info.txt')
	badg_setts = girvan_newman(badger_graph)
	post_gn_graph(badger_graph,badg_setts,8,'Doak - HW6 Badger Graph',graphspace,actual_badg)

if __name__ == "__main__":
	main()
