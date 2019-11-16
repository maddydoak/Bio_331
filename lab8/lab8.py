"""
Maddy Doak
Bio 331 Lab 8
Neighbor-Joining Algorithm by Saitou and Nei
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

def post_gn_graph(tree,title,gs_session):
	G = GSGraph()
	G.set_name(title)
	G.set_tags(['Lab8'])
	nodes = []
	for pair in tree.keys():
		if pair[0] not in nodes:
			nodes.append(pair[0])
		if pair[1] not in nodes:
			nodes.append(pair[1])
	for node in nodes:
		G.add_node(node,label=node)
		if '+' in node:
			G.add_node_style(node,
							 color = 'pink',
							 height = 30,
							 width = 30+20*(float(len(node))-1))
		else:
			G.add_node_style(node,
							 color = '#98edfa',
						 	 height = 30,
						 	 width = 30)
	for edge,branch_dist in tree.items():
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],width=float(branch_dist))
	G.set_data(data={'description': 'original OTUs are blue, all others are pink; edge thickness is proportional to branch length'})
	post(G, gs_session)

# Inputs: text file with pairwise distances between OTUs
# Return: distance matrix D
def read_dist_matrix(filename):
	D = {}
	matrix_rows = []
	with open(filename, 'rt') as file:
		for line in file:
			row = line.strip().split()
			matrix_rows.append(row)
	for node in matrix_rows[0]:
		D[node] = {}
		for row in matrix_rows[1:]:
			node2 = row[0]
			if node != node2:
				node_2_dist = row[matrix_rows[0].index(node)+1]
				D[node][node2] = float(node_2_dist)	# {node:{node2:distance}}
	return D

# Implementation of Neighbor-Joining Algorithm by Saitou and Nei
# Inputs: D = dictionary of distances between each pair of nodes: {node: [(node2, distance)]}
def neighbor_joining(D):
	T = {}
	while len(D.keys()) > 2:
		n = len(D.keys())

		# Making q
		q = {}
		banned = []
		for node in D.keys():						# node = i
			banned.append(node)
			for node2 in D[node].keys():			# node2 = j
				if node2 not in banned:				# if j hasn't previously been looked at
					ik_sum = 0
					for d in D[node].values():
						ik_sum += d
					jk_sum = 0
					for d in D[node2].values():
						jk_sum += d
					q[(node,node2)] = (n-2)*D[node][node2] - ik_sum - jk_sum
		
		# Finding nodes i,j with minimum q
		min_q = (None,100000)
		for nodes,q in q.items():
			if q < min_q[1]:
				min_q = (nodes,q)
		node_i = min_q[0][0]
		node_j = min_q[0][1]
		merge_node = node_i+"+"+node_j		# merged node is node "i+j"


		# Finding the branch distance between i and "i+j"; and j and "i+j"
		ik_dist = 0
		for dist in D[node_i].values():
			ik_dist += dist
		jk_dist = 0
		for dist in D[node_j].values():
			jk_dist += dist
		iu_dist = 0.5*D[node_i][node_j] + (1/(2*(n-2)))*(ik_dist-jk_dist)
		ju_dist = D[node_i][node_j] - iu_dist

		# Adding new edges and their branch distance to T
		T[(node_i,merge_node)] = iu_dist 	# T[(i,u)] = distance
		T[(node_j,merge_node)] = ju_dist 	# T[(j,u)] = distance

		# Add "i+j" to matrix D with calculated distance from it to all other nodes
		D[merge_node] = {}
		for node_k in D.keys():
			if node_k != merge_node and node_k != node_i and node_k != node_j:	# if k != u,i,j
				D[merge_node][node_k] = 0.5*(D[node_i][node_k]+D[node_j][node_k]-D[node_i][node_j])
				D[node_k][merge_node] = D[merge_node][node_k]
		
		# Delete nodes i and j (and all references to them) from matrix D
		nodes_to_delete = []
		for node in D.keys():
			for node2 in D[node]:
				if node2 == node_i or node2 == node_j:
					nodes_to_delete.append((node,node2))
		for node_pair in nodes_to_delete:
			del D[node_pair[0]][node_pair[1]]
		del D[node_i]
		del D[node_j]

	# Update T with the distance between the final two nodes
	if len(D.keys()) == 2:
		nodei = None
		nodej = None
		for node in D.keys():
			for node2 in D[node].keys():
				nodei = node
				nodej = node2
				break
		T[(nodei,nodej)] = D[nodei][nodej]
	
	return T

def main():
	class_matrix = read_dist_matrix('class-example.txt')
	class_tree = neighbor_joining(class_matrix)
	post_gn_graph(class_tree,"Doak, Lab8 - Class Example",graphspace)

	nj_ex_matrix = read_dist_matrix('NJ-example.txt')
	nj_tree = neighbor_joining(nj_ex_matrix)
	post_gn_graph(nj_tree,"Doak, Lab8 - NJ Example",graphspace)

	wiki_matrix = read_dist_matrix('wikipedia-example.txt')
	wiki_tree = neighbor_joining(wiki_matrix)
	post_gn_graph(wiki_tree,"Doak, Lab8 - Wiki Example",graphspace)

if __name__ == "__main__":
	main()