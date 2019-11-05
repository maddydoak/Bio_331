# Maddy Doak
# HW1 - 9/18/2019
# BIOL 331

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph
import math # Only used for a log calculation for getting the k value

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def read_matrix(fileName):
	file = open(fileName, 'rt')
	nodes = file.readline().split()
	adjMatrix = []
	for line in file:
		listTimes = line.split()
		adjMatrix.append(listTimes[1:])
	return adjMatrix, nodes[1:]

def check_symmetric(adjMatrix):
	check = "Yes!"
	nRows = len(adjMatrix)
	while nRows > 0 and check == "Yes!":
		for row in range(len(adjMatrix)):
			for col in range(len(adjMatrix)-(nRows-1),len(adjMatrix)):
				# Starts at column 1, then 2, then 3, etc, so no redundancy
				if adjMatrix[row][col] != adjMatrix[col][row]:
					check == "No :("
			nRows -= 1
	return check

def mat_to_edgelist(adjMatrix,nodes):
	edgelist = []
	nRows = len(adjMatrix)
	for row in range(len(adjMatrix)):
		for col in range(len(adjMatrix)-(nRows-1),len(adjMatrix)):
			# Same as check_symmetric, starts at column 1, then 2, then 3, so no duplicates
			if adjMatrix[row][col] != '0':
				edgelist.append([nodes[row],nodes[col],adjMatrix[row][col]])
				# Adds node 1, node 2, and the amount of time spent together as an edge
		nRows -= 1
	return edgelist

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def post_graph(gs_session,nodes,edges):
	G = GSGraph()
	G.set_name('Badger Graph - Doak')
	G.set_tags(['HW1'])
	badgerInfo = parse_info("BadgerInfo.txt")
	degreeOfNode = [0]*len(nodes)
	# Create a list of the same length as nodes
	# Then check node 1 in every edge, find its index in nodes,
	# and add one to its degree in the degreeOfNode list
	for edge in edges:
		indBadger = nodes.index(edge[0])
		degreeOfNode[indBadger] += 1
	nodeCount = 0 # For the index of each node as it iterates through, for finding
				  # the degree of that node
	for node in nodes:
		if badgerInfo[node][1] == 'N':
			infoTwo = 'Negative for TB'
		else:
			infoTwo = 'Positive for TB'
		G.add_node(node, label = node, k = 1, popup = badgerInfo[node][0]+'<br>'+
													  infoTwo+'<br>'+
													  'Social group '+badgerInfo[node][2])
		borderStyle = 'solid'	# Throws error if not defined early in the function
		if badgerInfo[node][0] == 'Male':
			shape = 'rectangle'
		elif badgerInfo[node][0] == 'Female':
			shape = 'ellipse'
		else:
			shape = 'triangle'
		if badgerInfo[node][1] == 'P': 		# positive for TB
			color = '#1786a3'
		elif badgerInfo[node][1] == 'N':	# negative for TB
			color = '#97cddb'
		else:								# unknown TB status
			color = 'white'
		if badgerInfo[node][2] in ['1','2','3']:
			borderColor = 'black'
			if badgerInfo[node][2] == '1':
				borderWidth = 2
			elif badgerInfo[node][2] == '2':
				borderWidth = 5
			elif badgerInfo[node][2] == '3':
				borderWidth = 9
		elif badgerInfo[node][2] in ['4','5','6']:
			borderColor = 'black'
			borderStyle = 'dashed'
			if badgerInfo[node][2] == '4':
				borderWidth = 2
			elif badgerInfo[node][2] == '5':
				borderWidth = 5
			elif badgerInfo[node][2] == '6':
				borderWidth = 9
		elif badgerInfo[node][2] in ['7','8']:
			borderColor = 'gray'
			if badgerInfo[node][2] == '7':
				borderWidth = 2
			elif badgerInfo[node][2] == '8':
				borderWidth = 5
		sizeMod = degreeOfNode[nodeCount]*3 # Degree * 3 to add to the size of each node
		G.add_node_style(node, height = 40+sizeMod, width = 60+sizeMod, shape = shape, color = color, 
						 border_color = borderColor, border_width = borderWidth, style = borderStyle)
		nodeCount += 1
	for edge in edges:
		G.add_edge(edge[0],edge[1],k=10/(math.log(int(edge[2]),10)+0.1),
				   popup='Number of seconds interacted: '+edge[2])
		G.add_edge_style(edge[0],edge[1],
						 width = int(edge[2])/10000 + 0.5)
	G.set_data(data = {'description': 'Males = rectangles, \
									   females = ellipses, \
									   unknown = triangles;<br>\
									   TB positive = dark blue, \
									   TB negative = light blue, \
									   unknown TB status = white;<br>\
									   social group 1-3 = black border (width 2,5,9);<br>\
									   social group 4-6 = black dashed border  (width 2,5,9);<br>\
									   social group 7-8 = gray border (width 2,5);<br>\
									   size of node indicates degree;<br>\
									   size of edge indicates amount of time spent together'})
	post(G, gs_session)

def parse_info(fileName):
	badgerDict = {}
	file = open(fileName, 'rt')
	for line in file:
		wordList = line.split()
		attributes = ()
		for word in wordList[1:]:
			attributes += (word,)
		badgerDict.update({wordList[0] : attributes})
	return badgerDict

def main():
	adjMatrix,nodes = read_matrix("BadgerMatrix.txt")
	print("Is the adjacency matrix symmetric? Checking...")
	print(check_symmetric(adjMatrix))
	edges = mat_to_edgelist(adjMatrix,nodes)
	post_graph(graphspace,nodes,edges)

if __name__ == '__main__':
  main()