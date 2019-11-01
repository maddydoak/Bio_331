"""
Maddy Doak
BIOL 331 - Lab 6 ("Network Motifs")
10/30/2019

The goal of the lab is to find network motifs in a directed graph and verify if 
they are statistically significant
"""

"""
Format of .txt files:
node1;node2
"""

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph
import random

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def get_graph_dict(filename):
	graph = {}
	with open(filename, 'rt') as file:
		for line in file:
			line = line.strip()
			nodes = line.split(';')
			if nodes[0] not in graph.keys():
				graph[nodes[0]] = [nodes[1]]
			else:
				if nodes[1] not in graph[nodes[0]]:
					graph[nodes[0]].append(nodes[1])
	return graph 										# {node: [outgoing edges]}

# Helper function for count_motifs
def walk_back(s, graph):								# s = starting node (string)
	motifs = {'FFL':[0,[]],
			  'MCL':[0,[]],
			  'SIM':[0,[]]}
	if len(graph[s]) > 1:
		SIM_count = 0
		SIM_edges = []
		FFL_edges = []
		for s_nei in graph[s]:
			if s_nei in graph.keys():
				for nei in graph[s_nei]:
					if nei in graph[s]:
						if s_nei != nei:
							edges = [[s,s_nei],[s_nei,nei],[s,nei]]
							motifs['FFL'][0] += 1
							motifs['FFL'][1].append(edges)
			temp_count = 0
			for node, outgoing in graph.items():
				if s_nei in outgoing:
					temp_count += 1
			if temp_count == 1:
				SIM_count += 1
				SIM_edges.append([s,s_nei])
		if SIM_count == 2:
			motifs['SIM'][0] += 1
			motifs['SIM'][1] = SIM_edges
	visited = [s]
	neighbors = []
	for s_nei in graph[s]:
		neighbors.append(s_nei)
	predecessor = [s]
	MCL_edges = []
	while len(neighbors) > 0:
		curr = neighbors.pop()
		visited.append(curr)
		if curr in graph.keys():
			predecessor.append(curr)
			for nei in graph[curr]:
				if nei not in visited:
					neighbors.append(nei)
				if nei == s:
					for i in range(len(predecessor[:-1])):
						MCL_edges.append([predecessor[i],predecessor[i+1]])
					motifs['MCL'][0] += 1
					for edge in MCL_edges:
						if edge not in motifs['MCL'][1]:
							motifs['MCL'][1].append(edge)
	return motifs

def count_motifs(graph):
	# FFL: origin node connects directly to a target and to an intermediate that connects to that target
	# MCL: A cycle
	# SIM: origin node connects to 2+ targets that have no other incoming neighbors 
	motifs_edges = {'FFL':[0,[]],
			  		'MCL':[0,[]],
			  		'SIM':[0,[]]}	
	for start_node, outgoing in graph.items():
		curr_motifs = walk_back(start_node,graph)
		for motif in motifs_edges.keys():
			motifs_edges[motif][0] += curr_motifs[motif][0]
			for edge in curr_motifs[motif][1]:
				if edge not in motifs_edges[motif][1]:
					motifs_edges[motif][1].append(edge)
	return motifs_edges

def rewire_graph(graph,r):	# r = number of rewirings (int)
	nodes = []
	for node, outgoing in graph.items():
		if node not in nodes:
			nodes.append(node)
		for out in outgoing:
			if out not in nodes:
				nodes.append(out)
	for i in range(r):
		node1 = random.choice(nodes)
		while node1 not in graph.keys():
			node1 = random.choice(nodes)
		node2 = random.choice(nodes)
		while node2 not in graph.keys():
			node2 = random.choice(nodes)
		while node2 == node1 or node2 in graph[node1] or node1 in graph[node2]:
			node2 = random.choice(nodes)
			while node2 not in graph.keys():
				node2 = random.choice(nodes)
		node1_out = random.choice(graph[node1])
		while node1_out in graph[node2]:
			node1_out = random.choice(graph[node1])
		node2_out = random.choice(graph[node2])
		while node2_out in graph[node1]:
			node2_out = random.choice(graph[node2])
		graph[node1].remove(node1_out)
		graph[node1].append(node2_out)
		graph[node2].remove(node2_out)
		graph[node2].append(node1_out)
	return graph

def get_p_value(graph,r,t):	# r = number of rewirings, t = number of iterations
	orig_motifs = count_motifs(graph)
	original = orig_motifs['FFL'][0]
	pval = 0
	for i in range(t):
		new_graph = rewire_graph(graph,r)
		new_motifs = count_motifs(new_graph)
		if new_motifs['FFL'][0] >= original:
			pval += 1
	return pval/t

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def post_motif_graph(graph,motifs,gs_session,name,pval):
	nodes = []
	edges = []
	for node, outgoing in graph.items():
		if node not in nodes:
			nodes.append(node)
		for out in outgoing:
			if out not in nodes:
				nodes.append(out)
			if [node,out] not in edges:
				edges.append([node,out])
	G = GSGraph()
	G.set_name(name)
	G.set_tags(['Lab6'])
	for node in nodes:
		G.add_node(node,label=node)
		G.add_node_style(node,
						 color = '#56b5bf', 
						 border_color = '#3b7b82',
						 border_width = 2,
						 height = 30,
						 width = 30)
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		if edge in motifs['FFL'][1][0]:
			G.add_edge_style(edge[0],edge[1],
							 directed = True,
							 color = 'red')
		else:
			G.add_edge_style(edge[0],edge[1],
						 directed = True)
	G.set_data(data={'description': 'Motif is FFL. Freq = '+str(motifs['FFL'][0])+', r = 10, t = 5, p-value = '+str(pval)})
	post(G, gs_session)

def main():
	file_name = 'example.txt'
	graph = get_graph_dict(file_name)
	dict_motifs = count_motifs(graph)
	#post_motif_graph(graph,dict_motifs,graphspace)

	big_file = 'yeast-GRN-full.txt'
	big_graph = get_graph_dict(big_file)
	big_dict_motifs = count_motifs(big_graph)
	pval = get_p_value(big_graph,10,5)
	post_motif_graph(big_graph,big_dict_motifs,graphspace,'Doak - FFL Motif Yeast Graph',pval)


if __name__ == "__main__":
	main()