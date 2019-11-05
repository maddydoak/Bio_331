# Maddy Doak
# HW2, BIOL 331

import matplotlib.pyplot as plt
import numpy as np
import math

DUMMY_PATH_LEN = 10000

# Inputs: file name (string)
# Outputs: edge list (list or set of tuples)
def read_edge_file(file_name):
	edges = []
	nodes_degs = {}								# keeping track of number of nodes, deg(node), neighbors(node)
	with open (file_name, 'rt') as file:
		for line in file:
			edges.append(line.split())
			for node in line.split():
				if node not in nodes_degs:
					nodes_degs.update({node : [1,[]]})		# node_name : [deg, [neighbors]]
				else:
					nodes_degs[node][0] += 1				# updating degree for known nodes
			nodes = [line.split()[0],line.split()[1]]
			node_neighbors = [nodes_degs[nodes[0]][1], nodes_degs[nodes[1]][1]]
			for i in range(len(nodes)):
				if nodes[i] not in node_neighbors[-(i+1)]:	# updating list of neighbors if not already included
					node_neighbors[-(i+1)].append(nodes[i])
	return edges, len(nodes_degs), len(edges), nodes_degs	# Other outputs are for report

# uses list of all node degrees to plot a histogram of 
# node degree distribution
# matplotlib code taken from Anna's example
def plot_deg_distrib(node_values,ds_name,count):
	list_node_deg = {}
	for node, values in node_values.items():	# pulls degrees from dictionary of nodes/degrees/neighbors
		if values[0] not in list_node_deg.keys():
			list_node_deg[values[0]] = 1
		else:
			list_node_deg[values[0]] += 1
	x = [key for key in list_node_deg.keys()]
	y = [val for val in list_node_deg.values()]
	fig_n = str(count)+'0'
	plt.figure(int(fig_n))
	plt.bar(x,y)
	plt.xlim([1,100])
	plt.xscale('log')
	plt.yscale('log', nonposy='clip') 			# from https://stackoverflow.com/questions/17952279/logarithmic-y-axis-bins-in-python
	plt.xlabel('log(degree (k))')
	plt.ylabel('log(count)')
	plt.title('Degree Distribution')
	plt.savefig('deg_distrib_'+ds_name+'.png')
	print('wrote to deg_distrib_'+ds_name+'.png')

# matplotlib code adapted from Anna's example
def plot_avg_AND_C(edge_list,node_values,ds_name):
	deg_avg_AND_C = {}
	for node, values in node_values.items():
		deg = int(values[0])
		if deg not in deg_avg_AND_C:		# Keeps a dictionary of degree k and its nodes' ANDs then C(v)s
			deg_avg_AND_C[deg] = [[],[]]	# { degree : [[node ANDs],[node C(v)s]]}
		sum_neighbor_deg = 0
		clust_coef = 0
		neighbors = values[1]
		if deg == 1:
			# Sets sum_neighbor_degree to degree of only neighbor
			sum_neighbor_deg = node_values[neighbors[0]][0]
		elif deg != 0:
			neigh_edges = 0						# Unique edges between neighbors of node
			for neighbor in neighbors:			# Adds unique edges among neighbors to list
				sum_neighbor_deg += node_values[neighbor][0]
				meta_neighbors = node_values[neighbor][1]
				for meta in meta_neighbors:
					if meta in neighbors:
						neigh_edges += 1
			clust_coef = neigh_edges / (deg*(deg-1))
		deg_avg_AND_C[deg][0].append(sum_neighbor_deg / len(neighbors))
		deg_avg_AND_C[deg][1].append(clust_coef)

	for deg, list_AND_C in deg_avg_AND_C.items():	# Averages ANDs and C(v)s for nodes by degree
		avg_AND = 0
		avg_C = 0
		for AND in list_AND_C[0]:
			avg_AND += AND
		for C in list_AND_C[1]:
			avg_C += C
		avg_AND = avg_AND / len(list_AND_C[0])
		avg_C = avg_C / len(list_AND_C[1])
		deg_avg_AND_C[deg] = [avg_AND,avg_C]
	deg_avg_AND_C = sorted(deg_avg_AND_C.items())

	# Plotting average AND by degree
	x1 = [deg_AND_C[0] for deg_AND_C in deg_avg_AND_C]
	y1 = [deg_AND_C[1][0] for deg_AND_C in deg_avg_AND_C]
	plt.figure(1)
	plt.scatter(x1, y1, label=ds_name)			# from https://pythonspot.com/matplotlib-bar-chart/
	plt.xlim([1,100])
	plt.legend()
	plt.xscale('linear')
	plt.yscale('linear')
	plt.xlabel('Degree of Node (K)')
	plt.ylabel('Average AND')
	plt.title('Average AND by Node Degree')
	plt.savefig('avg_AND_'+ds_name+'.png')
	print('wrote to avg_AND_'+ds_name+'.png')

	# Plotting average C(v) by degree
	x2 = [deg_AND_C[0] for deg_AND_C in deg_avg_AND_C]
	y2 = [deg_AND_C[1][1] for deg_AND_C in deg_avg_AND_C]
	plt.figure(2)
	plt.scatter(x2, y2, label=ds_name)							# from https://pythonspot.com/matplotlib-bar-chart/
	plt.xlim([1,100])
	plt.legend()
	plt.xscale('linear')
	plt.xlabel('Degree of Node (K)')
	plt.ylabel('Average C(v)')
	plt.title('Average C(v) by Node Degree')
	plt.savefig('avg_C_'+ds_name+'.png')
	print('wrote to avg_C_'+ds_name+'.png')


# From lab4 (written by myself, modified slightly for HW2)
# Returns dictionary of {neighbor (str) : len_shortest_path (int)}
def shortest_paths(node_values,banned_nodes,starting_node):
	distances = {}
	for node in node_values.keys():
		if node not in banned_nodes:
			distances[node] = DUMMY_PATH_LEN
	distances[starting_node] = 0
	visited_nodes = [starting_node]
	while len(visited_nodes) != 0:
		curr = visited_nodes.pop(0)
		if curr in node_values.keys():
			for neighbor in node_values[curr][1]:
				if neighbor in distances:
					if distances[neighbor] == DUMMY_PATH_LEN:
						distances[neighbor] = distances[curr] + 1
						visited_nodes.append(neighbor)
				else:
					distances[neighbor] = distances[curr] + 1
					visited_nodes.append(neighbor)
	return distances

def bfs_hist(node_values,ds_name):
	n_paths_len_k = {}							# Dict of {path_len : n_paths}
	banned_nodes = []
	for node in node_values.keys():
		distances = shortest_paths(node_values,banned_nodes,node)
		for neighbor, path_len in distances.items():
			if path_len != DUMMY_PATH_LEN and path_len != 0:
				if path_len not in n_paths_len_k.keys():
					n_paths_len_k[path_len] = 0
				n_paths_len_k[path_len] = n_paths_len_k[path_len] + 1
		banned_nodes.append(node)
	n_paths_len_k = sorted(n_paths_len_k.items())

	x = [paths[0] for paths in n_paths_len_k]
	y = [paths[1] for paths in n_paths_len_k]
	y_pos = np.arange(x[-1])
	plt.figure(3)
	plt.scatter(x, y, label=ds_name)							# from https://pythonspot.com/matplotlib-bar-chart/
	plt.xticks(x, x)
	plt.xscale('linear')
	plt.xlim([1,100])
	plt.legend()
	plt.xlabel('Length of Shortest Path Between Nodes')
	plt.ylabel('Number of Paths')
	plt.title('Distribution of Node Pairs with Each Possible Shortest Path Length')
	plt.savefig('path_len_'+ds_name+'.png')
	print('wrote to path_len_'+ds_name+'.png')

#############################################################################

def main():
	names = ['Fly','Yeast-LC','HIPPIE','Yeast-Y2H','Yeast-APMS']
	files = ['input-files/Fly_Unpublished.txt', 'input-files/Yeast_LC_Multiple.txt',
			'input-files/HIPPIE_Unweighted.txt', 'input-files/Yeast_Y2H_Union.txt',
			'input-files/Yeast_Combined_APMS.txt']
	data_dic = {}
	for i in range(len(names)):
		print('DATASET:',names[i])
		print('READING FILE:',files[i])
		# node_info = dict with entries: {node_name (str) : [deg (int), [neighbors] (list of str)]}
		edges, n_nodes, n_edges, node_info = read_edge_file(files[i])
		data_dic[names[i]] = [edges,n_nodes,n_edges,node_info]
	for dataset, info in data_dic.items():
		print(dataset)
		print(info[1],info[2])
	count = 1
	for dataset, info in data_dic.items():
		print('Plotting information for: '+dataset)
		plot_deg_distrib(info[3],dataset,count)
		count += 1
		plot_avg_AND_C(info[0],info[3],dataset)
	yeast = ['Yeast-LC','Yeast-Y2H','Yeast-APMS']
	yeast_files = ['input-files/Yeast_LC_Multiple.txt','input-files/Yeast_Y2H_Union.txt',
			'input-files/Yeast_Combined_APMS.txt']
	for i in range(len(yeast)):
		bfs_hist(data_dic[yeast[i]][3],yeast[i])

if __name__ == '__main__':
  main()