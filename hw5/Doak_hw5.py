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
	with open(filename, 'rt') as file:
		for line in file:
			edge = line.strip().split()
			for node in edge[:2]:
				if node not in nodes:
					nodes.add(node)
			edges.append(edge)
	return nodes, edges

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
		min_node = ['NA',100000]		# [node name, minimum distance value]
		for node in visited:
			if distances[node] < min_node[1]:
				min_node = [node,distances[node]]
		neighbors = []
		for edge in w_edges:
			if edge[0] == min_node[0] and edge[1] != predecessor[edge[0]]:
				neighbors.append([edge[1],edge[2]])
			elif edge[1] == min_node[0] and edge[0] != predecessor[edge[1]]:
				neighbors.append([edge[0],edge[2]])
		for neighbor in neighbors:
			if distances[neighbor[0]] == INITIAL_DISTANCE:
				visited.append(neighbor[0])
			temp_dist = float(distances[min_node[0]]) + float(neighbor[1])
			if temp_dist < distances[neighbor[0]]:
				distances[neighbor[0]] = temp_dist
				predecessor[neighbor[0]] = min_node[0]
		visited.remove(min_node[0])
	return distances, predecessor

# Recursive helper function for get_path
def get_path_helper(pi,u,path):
	if pi[path[0]] != 0:
		pred = pi[path[0]]
		path.insert(0,pred)
		if pred == u:
			return path
		else:
			return get_path_helper(pi,u,path)
	else:
		return path

# Gets shortest path from node u to node v
# given predecessors (pi)
# Returns list of nodes in that path
def get_path(pi,u,v):
	if pi[u] == v or pi[v] == u:
		return [u,v]
	else:
		return get_path_helper(pi,u,[v])

# Returns the edges in a Steiner tree using the metric closure and it's minimum spanning tree
def steiner_approx(nodes,w_edges,terminal_file):
	t = []	# terminals
	with open(terminal_file, 'rt') as file:
		for line in file:
			terminal = line.strip()
			if terminal not in t:
				t.append(terminal)
	metric_closure = []							# {(terminal1,terminal2) = length(shortest path)}
	for i in range(len(t)-1):
		dist,pi = dijkstra(nodes,w_edges,t[i])
		for j in range(i+1,len(t)):
			path = get_path(pi,t[i],t[j])
			path_weight = 0
			for k in range(len(path)-1):
				for edge in w_edges:
					if path[k] in edge and path[k+1] in edge:
						path_weight += float(edge[2])
			metric_closure.append([t[i],t[j],path_weight])
	MC_MST_edges = kruskal(metric_closure)		# List of edges in metric closure's MST
	tree = []
	for edge in MC_MST_edges:
		u = edge[0]
		v = edge[1]
		dist,pi = dijkstra(nodes,w_edges,u)
		p = get_path(pi,u,v)
		nodes_p_in_t = []
		for node in p:
			for edge in tree:
				if node in edge:
					nodes_p_in_t.append(node)
		if len(nodes_p_in_t) < 2:
			for i in range(len(p)-1):
				for edge in w_edges:
					if p[i] in edge and p[i+1] in edge:
						tree.append([edge[0],edge[1],edge[2]])
		else:
			p_i = nodes_p_in_t[0]
			p_j = nodes_p_in_t[-1]
			pi_path = get_path(pi,u,p_i)
			pj_path = get_path(pi,p_j,v)
			for path in [pi_path,pj_path]:	# Add pi_path and pj_path to T
				for i in range(len(path)-1):
					for edge in w_edges:
						if path[i] in edge and path[i+1] in edge:
							tree.append([edge[0],edge[1],edge[2]])
	return tree

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def graph_steiner(nodes,edges,steiner_edges,terminal_file,gs_session,name):
	G = GSGraph()
	G.set_name(name)
	G.set_tags(['HW5'])
	t = []	# terminals
	with open(terminal_file, 'rt') as file:
		for line in file:
			terminal = line.strip()
			if terminal not in t:
				t.append(terminal)
	for node in nodes:
		G.add_node(node,label=node)
		if node not in t:
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
		if edge in steiner_edges:
			G.add_edge_style(edge[0],edge[1],width=str(float(edge[2])+1),
							 color='red')
		else:
			G.add_edge_style(edge[0],edge[1],width=str(float(edge[2])+1))
	for edge in steiner_edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],width=str(float(edge[2])+1),
							 color='red')
	G.set_data(data={'description': 'Steiner tree edges in red, terminals in red'})
	post(G, gs_session)

def main():
	ex_graph = 'example-graph.txt'
	node_set,edge_list = get_nodes_edges(ex_graph)
	steiner_edges = steiner_approx(node_set,edge_list,'example-terminals.txt')
	graph_steiner(node_set,edge_list,steiner_edges,'example-terminals.txt',graphspace,'Doak HW5 Test Graph')

	big_graph = 'EGFR1-graph.txt'
	big_nodes,big_edges = get_nodes_edges(big_graph)
	steiner_big = steiner_approx(big_nodes,big_edges,'EGFR1-terminals.txt')
	graph_steiner(big_nodes,big_edges,steiner_big,'EGFR1-terminals.txt',graphspace,'Doak HW5 EGFR1 Graph')

if __name__ == "__main__":
	main()