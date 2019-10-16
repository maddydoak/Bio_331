"""
Maddy Doak
BIOL 331 - HW4 ("Minimum Spanning Tree")
main() function at bottom

Implementing Kruskal's Algorithm for finding the Minimum Spanning Tree (MST) of a graph
"""

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

# Input: list of edges (lists of 2 nodes with optional additional info) 
# Returns: list of nodes
def get_nodes(edges):
	nodes = []
	for edge in edges:
		if edge[0] not in nodes:
			nodes.append(edge[0])
		elif edge[1] not in nodes:
			nodes.append(edge[1])
	return nodes

def post_MST_graph(all_edges,MST_edges,name,list_tags,gs_session):
	nodes = get_nodes(all_edges)
	G = GSGraph()
	G.set_name(name)
	G.set_tags(list_tags)
	for node in nodes:
		G.add_node(node,label=node)
		G.add_node_style(node,
						 color = '#56b5bf',
						 height = 50,
						 width = 50)
	for edge in all_edges:
		G.add_edge(edge[0],edge[1],popup='weight: '+str(edge[2]))
		if edge in MST_edges:
			G.add_edge_style(edge[0],edge[1],width=edge[2],color='#6edfeb')
		else:
			G.add_edge_style(edge[0],edge[1],edge_style='dotted',width=edge[2],color='#fc8f81')
	G.set_data(data = {'description': 'Popups show the weights of the edges.\
										Edges included in the MST are colored\
										in blue, all excluded edges are in red.'})
	post(G, gs_session)

def post_EGFR_graph(MST_edges,name,list_tags,gs_session):
	nodes = get_nodes(MST_edges)
	G = GSGraph()
	G.set_name(name)
	G.set_tags(list_tags)
	for node in nodes:
		G.add_node(node,label=node)
		G.add_node_style(node,
						 color = '#56b5bf',
						 height = 50,
						 width = 50)
	for edge in MST_edges:
		G.add_edge(edge[0],edge[1],popup='weight: '+str(edge[2]))
		G.add_edge_style(edge[0],edge[1],width=edge[2],color='#6edfeb')
	G.set_data(data = {'description': 'Popups show the weights of the edges.\
										Only MST edges are included.'})
	post(G, gs_session)

# Reads file of weighted, undirected edges with three-columns (node1, node2, and weight)
def read_edge_file(file_name):
	edges = []
	with open (file_name, 'rt') as file:
		for line in file:
			edges.append(line.split())
	return edges							# edges: [[node1,node2,weight]]

# Implementing Kruskal's Algorithm for finding the Minimum Spanning Tree (MST) of a graph
# Inputs: weighted, undirected graph (list of edges, each edge is [node1,node2,weight]) 
# Returns: edges in the MST (list of edges)
def kruskal(edges):
	forest = {}								# {(nodes in tree) : [edges in tree]}
	nodes = get_nodes(edges)
	for node in nodes:
		forest[(node,)] = []
	for edge in edges:
		edge[2] = int(float(edge[2]))
	edges = sorted(edges, key=lambda edge: edge[2])
	while len(edges) > 0:
		curr = edges.pop(0)					# Pops out the edge with minimum weight from edges
		is_tree = True
		tree1 = -1
		tree2 = -1
		while is_tree and (tree1 == -1 or tree2 == -1):
			for nodes in forest.keys():
				if curr[0] in nodes:
					if curr[1] not in nodes:
						tree1 = nodes
					else:
						is_tree = False
				else:
					if curr[1] in nodes:
						tree2 = nodes
		if tree1 != tree2 and is_tree:
			# combine trees
			combined_tree = []
			for tree in [tree1,tree2]:
				for node in tree:
					combined_tree.append(node)
			combined_tree = tuple(combined_tree)
			forest[combined_tree] = forest[tree1]	# Combined tree starts with list of edges of tree1
			for edge in forest[tree2]:				# Then each edge from tree2 is appended
				forest[combined_tree].append(edge)
			for tree in [tree1,tree2]:
				del forest[tree]
			# add edge to tree
			forest[combined_tree].append(curr)
	MST_edges = []
	for value in forest.values():
		for edge in value:
			MST_edges.append(edge)
	return MST_edges

####################################################################################################################

def main():
	"""
	sample_edges = read_edge_file('weighted-graph.txt')
	MST_sample_edges = kruskal(sample_edges)
	post_MST_graph(sample_edges,MST_sample_edges,'HW4 - Weighted Sample Graph',['HW4'],graphspace)
	"""

	EGFR_edges = read_edge_file('weighted-EGFR.txt')
	MST_EGFR_edges = kruskal(EGFR_edges)
	post_EGFR_graph(MST_sample_edges,'HW4 - EGFR MST',['HW4'],graphspace)

	# Started around 8:46pm

if __name__ == "__main__":
	main()