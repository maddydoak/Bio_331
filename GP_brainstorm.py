"""
Maddy Doak
Brainstorm for Group Project
11/27/2019
"""

"""
IDEAS:
	Randomly pick unlabeled
	Node2vec
		Take graph, map each node to a vector; take node vector similarity
	Degree
		Find nodes with similar degrees, 
	Shared neighbors
		Find nodes with lots of shared neighbors
	Node centralities
		Look at centralities of positives for NMII, find other nodes with similar centralities
	Pattern of pos/neg neighbors
		Define a pattern of nodes with a number of pos/neg neighbors (use gene ontology)
	Steiner tree
		Terminals are pos nodes; then random walk to rank candidates
	Dijikstra-all
		Dij-all (all tied shortest paths) from Fog to Sqh (source --> target)
		Next step - remove "critical" nodes?
		Or find "k" shortest paths (1st, 2nd, 3rd shortest paths - etc) 
	Yen's K shortest paths - KSP
	Model Steiner Tree --> look at neighbors of tree
	Look at subnetwork of predictions, rank by number of methods that predict it
	Look at the intersection of all candidate predictions from different methods
	Shortest paths among all pos nodes (pattern(s)?)

	(Assumptions: these are regulators - that regulators are "close" to each other)
"""

def main():
	print("GP brainstorm")

if __name__ == "__main__":
	main()