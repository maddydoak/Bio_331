# Maddy Doak
# Lab 3, BIOL 331
# 9/18/2019

# Group 3

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph
import random

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def isSymmetric(node1,node2,edgeList):
	if node1 == node2 or [node1,node2] in edgeList or [node2,node1] in edgeList:
		return True
	else:
		return False

def makeERGraph(n,m):
	nodes = []
	for num in range(n):
		nodes.append(str(num))
	edges = []
	for num in range(m):
		node1 = str(random.randrange(n))
		node2 = str(random.randrange(n))
		while isSymmetric(node1,node2,edges):
			node1 = str(random.randrange(n))
			node2 = str(random.randrange(n))
		edges.append([node1,node2])
	return nodes,edges

def genBAGraph(n,m,t):
	nodes = []
	for num in range(n):
		nodes.append(str(num))
	edges = []
	for i in range(n-1):
		edges.append([nodes[i],nodes[i+1]])
	for i in range(t):
		newNode = str(int(nodes[-1])+1)
		nodes.append(newNode)
		for j in range(m):
			newConnection = edges[random.randrange(len(edges))][random.choice([0,1])] 
			while isSymmetric(newNode,newConnection,edges):
				newConnection = edges[random.randrange(len(edges))][random.choice([0,1])] 
			edges.append([newNode,newConnection])
	return nodes,edges

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def postGraph(nodes,edges,name,gs_session):
	G = GSGraph()
	G.set_name(name)
	G.set_tags(['Lab3'])
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
	post(G,gs_session)

#############################################################

def main():
	erNodes, erEdges = makeERGraph(25,100)
	baNodes, baEdges = genBAGraph(3,2,50)
	postGraph(erNodes,erEdges,'ER Graph',graphspace)
	postGraph(baNodes,baEdges,'BA Graph',graphspace)

if __name__ == '__main__':
	main()