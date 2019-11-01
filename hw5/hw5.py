"""
Maddy Doak
BIOL 331 - HW5 ("Minimum Steiner Tree Approximation Algorithm")
11/4/2019

Blurb
"""
from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph
from hw4 import kruskal											# HW4 Min Spanning Tree algorithm

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

# Takes input graph (.txt) in format node1, node2, weight
# Returns set of nodes, list of edges with weights
def get_nodes_edges(filename):
	nodes = set()
	edges = []
	line_count = 0
	with open(filename, 'rt') as file:
		for line in file:
			if line_count == 0:
				start_node = line.strip().split()[0]
			line_count += 1
			edge = line.strip().split()
			for node in edge[:2]:
				if node not in nodes:
					nodes.add(node)
			edges.append(edge)
	return nodes, edges, start_node

# Take set of nodes, list of edges with weights, starting node s
# Returns dict of distances with (node,dist) pairs and a dict of predecessors with
# (node, predecessor_node) pairs
INITIAL_DISTANCE = 10000
def dijkstra(nodes,w_edges,s):
	distances = {}
	predecessor = {}
	for node in nodes:
		distances[node] = INITIAL_DISTANCE
		predecessor[node] = 0
	distances[s] = 0
	visited = [s]
	while len(visited) != 0:
		min_node = ['NA',100000]		# node with minimum distance value
		for node in visited:
			if distances[node] < min_node[1]:
				min_node = [node,distances[node]]
		visited.remove(min_node[0])
		neighbors = []
		for edge in w_edges:
			if edge[0] == min_node and edge[1] != predecessor[edge[0]]:
				neighbors.append([edge[1],edge[2]])
			elif edge[1] == min_node and edge[0] != predecessor[edge[1]]:
				neighbors.append([edge[0],edge[2]])
		for neighbor in neighbors:
			if distances[neighbor[0]] == INITIAL_DISTANCE:
				visited.append(neighbor)
			temp_dist = distances[min_node] + neighbor[1]
			if temp_dist < distances[neighbor]:
				distances[neighbor] = temp_dist
				predecessor[neighbor] = min_node
	return distances, predecessor

# Recursive helper function for get_path
def get_path_helper(pi,u,path):
	pred = pi[path[0]]
	path_with_pred = path.insert(0,pred)
	if pred == u:
		return path_with_pred
	else:
		return get_path_helper(pi,u,path_with_pred)

# Gets shortest path from node u to node v
# given predecessors (pi)
# Returns list of nodes in that path
def get_path(pi,u,v):
	if pi[u] == v or pi[v] == u:
		return [u,v]
	else:
		return get_path_helper(pi,u,[v])

def steiner_approx(nodes,w_edges,pi,terminal_file):
	t = set()	# terminals
	with open(terminal_file, 'rt') as file:
		for line in file:
			terminal = line.strip()
			if terminal not in t:
				t.add(terminal)
	metric_closure = {}
	MC_MST_edges = kruskal(w_edges)			# CHECK THIS WORKS
	tree = {}
	for edge in met_clos_MST_edges:
		u = edge[0]
		v = edge[1]
		p = get_path(pi,u,v)
		nodes_p_in_t = []
		for node in p:
			if node in tree.keys():
				nodes_p_in_t.append(node)
		if len(nodes_p_in_t) < 2:
			for node in p:
				if node not in tree.keys():
					tree[node] = 0
		else:
			pi = p[0]
			pj = p[-1]
			pi_path = get_path(pi,u,pi)
			pj_path = get_path(pi,pj,v)
			# Add pi_path and pj_path to T
	return tree

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def graph(nodes,edges,gs_session,name):
	G = GSGraph()
	G.set_name(name)
	G.set_tags(['Lab6'])
	for node in nodes:
		G.add_node(node,label=node)
		if node[0] == 'u':
			G.add_node_style(node,
							 color = '#56b5bf', 
							 border_color = '#3b7b82',
							 border_width = 2,
							 height = 30,
							 width = 30)
		else:
			G.add_node_style(node,
							 color = '#fca59f', 
							 border_color = '#f78b83',
							 border_width = 2,
							 height = 30,
							 width = 30)
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],width=edge[2])
	G.set_data(data={'description': 'description here'})
	post(G, gs_session)

def main():
	ex_graph = 'example-graph.txt'
	node_set,edge_list,start_node = get_nodes_edges(ex_graph)
	dists,pi = dijkstra(node_set,edge_list,start_node)
	steiner_approx = (node_set,edge_list,pi,'example-terminals.txt')
	#graph(node_set, edge_list, graphspace, 'Doak HW5 Test Graph')

if __name__ == "__main__":
	main()