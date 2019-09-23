# Maddy Doak
# HW2, BIOL 331

from graphspace_python.api.client import GraphSpace 
from graphspace_python.graphs.classes.gsgraph import GSGraph
import matplotlib.pyplot as plt
import math

graphspace = GraphSpace('maddydoak@gmail.com','mygraphspacepw')

def post(G, gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

# Inputs: file name (string)
# Outputs: edge list (list or set of tuples)
#		   (REPORT: needs # nodes, # edges, AND())
def readEdgeFile(fileName):

# Anna's code from HW2 instructions
def plotSomeNumbers():
	fig = plt.figure(figsize=(6.5,4)) 		# make a 6.5" wide by 4" tall figure.
	x = list(range(1,10)) 					# make a list of [1,2,3,...9]
	y = [xval/1.5+2 for xval in x] 			# y is a function of x
	logx = [math.log(a) for a in x] 		# compute the log of x
	logy = [math.log(b) for b in y] 		# compute the log of y

	plt.subplot(1,2,1) 						# in a 1-by-2 grid, 1st subplot
	plt.plot(x,y,'o-r') 					# plot a red line with circles
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Some Numbers')

	plt.subplot(1,2,2) 						# in a 1-by-2 grid, 2nd subplot
	plt.plot(logx,logy,'s-b') 				# plot a blue line with squares
	plt.xlabel('log x')
	plt.ylabel('log y')
	plt.title('Some Numbers (log)')

	plt.tight_layout() 						# make the labels "snap" to the grid.
											# this may emit a warning, which is OK

	plt.savefig('numbers.png') 				# save figure as PNG
	print('wrote to numbers.png')
	return

def main():
	names = ['Fly','Yeast-LC','HIPPIE','Yeast-Y2H','Yeast-APMS']
	files = ['input-files/Fly_Unpublished.txt', 'input-files/Yeast_LC_Multiple.txt',
			'input-files/HIPPIE_Unweighted.txt', 'input-files/Yeast_Y2H_Union.txt',
			'input-files/Yeast_Combined_APMS.txt']
	for i in range(len(names)):
		print('DATASET:',names[i])
		print('READING FILE:',files[i])
		# use example network first!


################################################################

if __name__ == '__main__':
  main()