## Lab1 -- Warm Up
## September 4 2019
## Requirements: Python3
## By Maddy Doak

def main():
  """
  Main function. There are no inputs 
  (nodes and matrix are specified within this function)
  """

  ## List of node names
  nodes = ['A','B','C','D','E','F']
  ## Adjacency matrix (assume this represents an undirected graph & has same order as nodes)
  adj_mat = [[0,1,0,1,0,1],[1,0,0,0,1,1],[0,0,0,1,0,0],[1,0,1,0,0,0],[0,1,0,0,1,1],[1,1,0,0,1,0]]
  #adj_mat_challenge = [[0,1,0,1,0,1],[1,0,0,0,1,1],[0,0,0,1,0,0],[1,0,1,0,0,0],[0,1,0,0,1,1],[1,1,0,0,1,0]]

  ## function calls
  print('\nAdjacency Matrix:')
  print_mat(nodes,adj_mat)
  print("Number of Edges from Adj. Matrix:",num_edges_from_mat(adj_mat))

  print('\nMaking Adjacency List...')
  adj_list = mat_to_list(nodes,adj_mat)

  print('\nAdjacency List:')
  print_list(adj_list)
  print("Number of Edges from Adj. List:",num_edges_from_list(adj_list))

  print('\nDone!')
  return # done with main function

def print_mat(nodes,adj_mat):
  """
  Prints the adjacency matrix
  Inputs: nodes (list of strings) and adj_mat (list of lists)
  Returns: nothing
  """
  topLine = "  "
  for node in nodes:
    topLine += node
  print(topLine)
  nNodes = len(nodes)     # Used to keep track of the number of rows/cols left
  while nNodes > 0:       # Stops the loop when hits the end of the list of nodes
    line = ""
    line += nodes[-nNodes] + " "
    for num in adj_mat[-nNodes]:
      line += str(num)
    print(line)
    nNodes -= 1
  return

def num_edges_from_mat(adj_mat):
  """
  Counts the number of edges from the adjacency matrix
  Inputs: adj_mat (list of lists)
  Returns: the number of edges (int)
  """
  nEdges = 0
  colsToCount = 0         # Starts counting one to the right as you go down each row
  for row in adj_mat:
    for col in row[colsToCount:]:
      if col == 1:
        nEdges += 1
    colsToCount += 1
  return nEdges

def mat_to_list(nodes,adj_mat):
  """
  Converts the adjacency matrix to an adjacency list
  Inputs: nodes (list of strings) and adj_mat (list of lists)
  Returns: adjacency list (dictionary of (node,neighbor list) pairs).
  """
  connections = []
  for row in adj_mat:
    newRow = []           # Creates a list with letters instead of 0/1
    for colIndex in range(len(row)):
      if row[colIndex] == 1:
        newRow.append(nodes[colIndex])
    connections.append(newRow)
  adj_list = {
  }
  for nodeIndex in range(len(nodes)):
    adj_list.update({nodes[nodeIndex] : connections[nodeIndex]})
    """ 
    From https://thispointer.com/python-how-to-add-append-key-value-pairs-in-dictionary-using-dict-update/
    """
  return adj_list

def print_list(adj_list):
  """
  Prints the adjacency list
  Inputs: adj_list (dictionary of (node,list) pairs)
  Returns: nothing
  """
  for key,value in adj_list.items():    
    """
    From https://stackoverflow.com/questions/5904969/how-to-print-a-dictionarys-key
    """
    line = key + ": "
    for i in value:
      line += i + " "
    print(line)
  return 

def num_edges_from_list(adj_list):
  """
  Counts the number of edges from the adjacency list
  Inputs: adj_list (dictionary of (node,list) pairs)
  Returns: the number of edges (int)
  """
  nEdges = 0
  for key,value in adj_list.items():
    for item in value:
      if item == key:
        nEdges += 1          # To fix the problem of self-connections
    nEdges += len(value)
  return int(nEdges/2)

"""
 Leave this is at the bottom of the file. Once all functions are loaded, then 
 main() is called UNLESS you are importing this file into another script. 
 See https://docs.python.org/3/library/__main__.html
"""
if __name__ == '__main__':
  main()