# Maddy Doak
# HW3 - Bio 331
# Due 10/9/2019

import random
import math
from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def rgb_to_hex(red,green,blue):
	maxHexValue = 255
	r = int(red*maxHexValue)
	g = int(green*maxHexValue)
	b = int(blue*maxHexValue)
	RR = format(r,'02x')
	GG = format(g,'02x')
	BB = format(b,'02x')
	return '#'+RR+GG+BB

def post_graph(gs_session,node_dict,node_visits):
	v_max = 0
	v_min = 10000
	for node,n_visits in node_visits.items():
		node_visits[node] = [n_visits,n_visits]
	for n_visits in node_visits.values():
		n_visits[1] += 1
		n_visits[1] = math.log(n_visits[1])
		if n_visits[1] > v_max:
			v_max = n_visits[1]
		elif n_visits[1] < v_min:
			v_min = n_visits[1]
	for n_visits in node_visits.values():
		n_visits[1] = (n_visits[1] - v_min) / (v_max - v_min)
	G = GSGraph()
	G.set_name('Doak - HW3 EGFR Graph')
	G.set_tags(['HW3'])
	for node in node_dict.keys():
		G.add_node(node,label=node,popup='Count = '+str(node_visits[node][0])
										 +', normalized score = '+str(node_visits[node][1]))
		G.add_node_style(node,
						 color = rgb_to_hex(0.8*(1-node_visits[node][1]),1-node_visits[node][1],1),
						 height = 60,
						 width = 60)
	for node,nei in node_dict.items():
		for outgoing in nei[1]:
			if node not in node_dict[outgoing][1]:
				G.add_edge(node,outgoing,directed=True)
				G.add_edge_style(node,outgoing,
						 	 	 directed=True)
			else:
				node_dict[outgoing][1].remove(node)
				G.add_edge(node,outgoing)
	post(G, gs_session)

####################################################################
# Generating the data for the graph

# Inputs: filename for .txt file with set of directed edges
# Returns: dictionary of nodes and their incoming/outgoing neighbors
def parse_edges(filename):
	graph = {}
	with open (filename, 'rt') as file:
		for line in file:
			node1 = line.split()[0]
			node2 = line.split()[1]
			if node1 not in graph.keys():
				graph[node1] = [[],[node2]]
			else:
				if node2 not in graph[node1][1]:
					graph[node1][1].append(node2)
			if node2 not in graph.keys():
				graph[node2] = [[node1],[]]
			else:
				if node1 not in graph[node2][0]:
					graph[node2][0].append(node1)
	return graph									# dictionary, {node : [[incoming1,incoming2],[outgoing1,outgoing2]]}

# Inputs: graph (G) and source node (s)
# Returns: graph (G2) with added edges for every node v with no outgoing edges
def add_edges(G,s):
	G2 = {}
	for node,value in G.items():
		G2[node] = value
		if len(value[1]) == 0:						# if no outgoing neighbors:
			value[1].append(value[0][0])			# select first incoming neighbor to be an outgoing neighbor
	return G2

# Inputs: a representation of the graph (G), source node (s), n timesteps (T), teleportation probability (q)
# Returns: n nodes * T table with probabilities at each time step t (cols) for every node v (rows)
def rw_probs(G,s,T,q):
	nodes = []										# nodes = [A,B,C] then
	for node in G.keys():							# T = 2 --> [[A0,A1,A2],[B0,B1,B2],[C0,C1,C2]]
		nodes.append(node)
	prob_table = []
	for i in range(len(nodes)):
		row = []
		if i == 0:
			row.append(1)
		else:
			row.append(0)
		for j in range(T):
			row.append(0)
		prob_table.append(row)
	for i in range(len(nodes)):
		for j in range(1,T+1):
			summation = 0
			if len(G[nodes[i]][0]) > 0:
				for incoming in G[nodes[i]][0]:
					if j > 0:
						summation += (prob_table[nodes.index(incoming)][j-1]) / (len(G[incoming][1]))
			p = q*summation + (1-q)*prob_table[i][0]
			if p == 0.0:
				p = int(p)
			prob_table[i][j] = p 					# Probability of visiting node v at timestep t
	return prob_table

# Inputs: graph (G), source node (s), num time steps (T), teleportation probability (q)
# Returns: dict of {node (v) : num times walker visits v (int)}
def rw_simulate(G,s,T,q):
	num_visits = {}
	for node in G.keys():
		num_visits[node] = 0
	curr_node = s
	for t in range(T):
		choice = random.random()
		if choice > 0.85:
			curr_node = s
		else:
			curr_node = random.choice(G[curr_node][1])
		num_visits[curr_node] += 1
	return num_visits

def rw_error(prob_table,num_visits):
	nodes = []
	for node in num_visits.keys():
		nodes.append(node)
	# Normalizing numbers of visits to between 0,1
	v_max = 0
	v_min = 10000
	for n_visits in num_visits.values():
		if n_visits > v_max:
			v_max = n_visits
		elif n_visits < v_min:
			v_min = n_visits
	for n_visits in num_visits.values():
		n_visits = (n_visits - v_min) / (v_max - v_min)
	error = 0
	for i in range(len(nodes)):
		error += (num_visits[nodes[i]] - (prob_table[i][-1]))**2
	return error

#######################################################################

def main():
	T_table = 1000
	T_walk = 100000
	probability = 0.95
	EGFR_source = 'EGF'

	init_EGFR_graph = parse_edges('EGFR1-edges.txt')
	EGFR_graph = add_edges(init_EGFR_graph,EGFR_source)
	EGFR_table = rw_probs(EGFR_graph,EGFR_source,T_table,probability)
	EGFR_visits = rw_simulate(EGFR_graph,EGFR_source,T_walk,probability)
	error = rw_error(EGFR_table,EGFR_visits)

	print('Squared error for EGFR: '+str(error**2))
	post_graph(graphspace,EGFR_graph,EGFR_visits)

if __name__ == '__main__':
	main()