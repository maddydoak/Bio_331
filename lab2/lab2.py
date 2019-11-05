# Maddy Doak
# Lab 2 - 9/11/2019
# BIOL 331

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def post_test_graph(gs_session):
	nodes = ['A','B','C','D','E']
	edges = [['A','B'], ['C','A'],
			 ['A','D'],	['D','E'],
			 ['B','E'], ['C','D'],
			 ['B','D']]
	G = GSGraph()
	G.set_name('Test Graph - Doak')
	G.set_tags(['Lab 2'])
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
		G.add_edge_style(edge[0],edge[1],
						 edge_style = 'dotted')
	post(G, gs_session)


def post_dolphin_network(gs_session):
	nodes = []
	males = []
	females = []
	unknown = []
	for name in [['males.txt',males],
				 ['females.txt',females],
				 ['unknown-sex.txt',unknown]]:
		with open(name[0], 'rt') as file:  		   	   # From https://www.computerhope.com/issues/ch001721.htm
			for line in file:
				nodes.append(line.split())
				name[1].append(line.split())
	sideF = []
	upsideDown = []
	for name in [['side-floppers.txt',sideF],
				 ['upside-down-lobtailers.txt',upsideDown]]:
		with open(name[0], 'rt') as file:
			for line in file:
				name[1].append(line.split())
	edges = []
	with open ('dolphin_edgelist.txt', 'rt') as file:
		for line in file:
			edges.append(line.split())			   	   # From https://www.geeksforgeeks.org/python-string-split/
	G = GSGraph()
	G.set_name('Dolphin Network - Doak')
	G.set_tags(['Lab 2'])
	for node_list in nodes:
		for node in node_list:
			G.add_node(node,label=node)
			G.add_node_style(node, 
							 border_color = '#3b7b82',
							 border_width = 2)
	for listNames in [[males,'rectangle'],
					  [females,'ellipse'],
					  [unknown,'triangle']]:
		for item in listNames[0]:
			for name in item:
				G.add_node_style(name,
						 		 shape = listNames[1])
	for listNames in [[sideF,'#52989c'],
					  [upsideDown,'#82eff5']]:
		for item in listNames[0]:
			for name in item:
				G.add_node_style(name,
						 		 color = listNames[1])
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edge[0],edge[1],
						 edge_style = 'dotted')
	G.set_data(data = {'description': 'males=rectangles; females=ellipses; unknown=triangles'})
	post(G, gs_session)


def main():
	post_test_graph(graphspace)
	post_dolphin_network(graphspace)

if __name__ == '__main__':
  main()