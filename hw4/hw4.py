"""
Maddy Doak
BIOL 331 - HW4 ("Minimum Spanning Tree")
main() function at bottom

Implementing Kruskal's Algorithm for finding the Minimum Spanning Tree (MST) of a graph
Pushes graph with colored edges of MST to GraphSpace under specified user
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

# Helper function that returns list of nodes from list of edges
# Input: list of edges (lists of 2 nodes with optional additional info) 
# Returns: list of nodes
def get_nodes(edges):
	nodes = []
	for edge in edges:
		if edge[0] not in nodes:
			nodes.append(edge[0])
		if edge[1] not in nodes:
			nodes.append(edge[1])
	return nodes

# Posts a graph with different edge styles depending on whether the edge
# is in the MST or not
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
			G.add_edge_style(edge[0],edge[1],width=6+edge[2],color='#6edfeb')
		else:
			G.add_edge_style(edge[0],edge[1],edge_style='dotted',width=6+edge[2],color='#fc8f81')
	G.set_data(data = {'description': 'Popups show the weights of the edges.\
										Edges included in the MST are colored\
										in blue, all excluded edges are in red.'})
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
	nodes = get_nodes(edges)				# Helper function that returns a list of nodes (strings)
	for node in nodes:
		forest[(node,)] = []				# Each tree is a tuple (a key in the dict) where the forest starts 
	for edge in edges:						# with each node as its own tree with an empty list of edges
		edge[2] = float(edge[2])			# Converts edge weights to floats for sorting
	edges = sorted(edges, key=lambda edge: edge[2])		# Sorts edges by ascending weight
	while len(edges) > 0:
		curr = edges.pop(0)					# Pops out the edge with minimum weight ('curr') from edges
		is_tree = True						# Condition that breaks while loop if edge would create a cycle
		tree1 = -1							# The tree containing node1 from the edge
		tree2 = -1							# The tree containing node2 from the edge
		while is_tree and (tree1 == -1 or tree2 == -1):
			for nodes in forest.keys():		# While no cycles would be created and both nodes have not been
				if curr[0] in nodes:		# found in the forest...
					if curr[1] not in nodes:
						tree1 = nodes
					else:
						is_tree = False
				else:
					if curr[1] in nodes:
						tree2 = nodes 		# ...find where the nodes exist in the forest
		if tree1 != tree2 and is_tree:		# If the nodes are in separate trees and there are no cycles:
			# combine trees
			combined_tree = []
			for tree in [tree1,tree2]:
				for node in tree:
					combined_tree.append(node)
			combined_tree = tuple(combined_tree)
			forest[combined_tree] = forest[tree1]	# Combined tree starts with list of edges from tree1
			for edge in forest[tree2]:				# Then each edge from tree2 is appended
				forest[combined_tree].append(edge)
			for tree in [tree1,tree2]:				# After combining tree1/tree2, the originals are deleted 
				del forest[tree]
			# add edge to tree
			forest[combined_tree].append(curr)
	MST_edges = []									# List of edges from the forest (MST)
	for value in forest.values():
		for edge in value:
			MST_edges.append(edge)
	return MST_edges

####################################################################################################################

def main():

	# Reads in small example graph, finds the edges in the MST, and posts to Graphspace with
	# edges colored according to whether they exist in the MST
	sample_edges = read_edge_file('weighted-graph.txt')
	MST_sample_edges = kruskal(sample_edges)
	post_MST_graph(sample_edges,MST_sample_edges,'HW4 - Weighted Sample Graph (TEST)',['HW4'],graphspace)

	# Reads in full EGFR graph, finds the edges in the MST, and posts to Graphspace with
	# edges colored according to whether they exist in the MST
	EGFR_edges = read_edge_file('weighted-EGFR.txt')
	MST_EGFR_edges = kruskal(EGFR_edges)
	post_MST_graph(EGFR_edges,MST_EGFR_edges,'HW4 - EGFR MST',['HW4'],graphspace)

if __name__ == "__main__":
	main()