# Maddy Doak
# Lab 4 - BIOL 331
# 9/25/2019
# "Shortest Paths"

import Lab4_utils
from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def make_graph(nodes,edges,distances,session):
	G = GSGraph()
	G.set_name('Shortest Paths Graph - Doak')
	G.set_tags(['Lab4'])
	for node in nodes:
		G.add_node(node,label=node)
		G.add_node_style(node, 
			color = rgb_to_hex(0.8-distances[node]/10,1-distances[node]/10,1),
			width = 40,
			height = 40)
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],
						 edge_style = 'dotted')
	G.set_data(data = {'description': '(description here)'})
	post(G, session)

def shortest_paths(nodes,adj_list,starting_node):
	distances = {}
	for node in nodes:
		distances[node] = 10000
	distances[starting_node] = 0
	visited_nodes = [starting_node]
	while len(visited_nodes) != 0:
		curr = visited_nodes.pop(0)
		for neighbor in adj_list[curr]:
			if distances[neighbor] == 10000:
				distances[neighbor] = distances[curr] + 1
				visited_nodes.append(neighbor)
	return distances

def rgb_to_hex(red,green,blue):
	maxHexValue = 255
	r = int(red*maxHexValue)
	g = int(green*maxHexValue)
	b = int(blue*maxHexValue)
	RR = format(r,'02x')
	GG = format(g,'02x')
	BB = format(b,'02x')
	return '#'+RR+GG+BB


def main():
	node_list, edge_list, adj_list, adj_mat = Lab4_utils.get_graph('lab')
	distances = shortest_paths(node_list,adj_list,node_list[0])
	print('\nDistances:')
	for key,value in distances.items():
		print(str(key)+' = '+str(value))
	make_graph(node_list,edge_list,distances,graphspace)

if __name__ == '__main__':
	main()