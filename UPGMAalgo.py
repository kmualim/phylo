import itertools

class Node(object): 

	def __init__(self, dist, left=None, right=None): 
		self.dist = dist
		self.left = left 
		self.right = right 
		self.nn = None
		self.dist_to_NN = None

	# prints all the nodes 
	def __str__(self): 
		#list = []
		return ",".join(str(node) for node in self.__iter__())

	# iterate through all the taxa in Node and children 
	def __iter__(self): 
		for child in ([self.left,self.right]): 
			if child is None: 
				pass 
			# get item in child
			elif isinstance(child, Node): 
				for item in child: 
					yield item 
			# get child
			else: 
				yield child

	def get_len(self, node):
		return len(str(node))

	# gets distances between two nodes 
	def get_dist(self, node):
		# iterate through pairs again and again
		pairs = itertools.product(self, node) 
		distance = sum(self.dist(*pair) for pair in pairs)
		return distance/(self.get_len(self) + self.get_len(node))

	# updates distances 
	def update_distances(self, clusters):
		self.nn = None
		self.dist_to_NN = None 
		for node in clusters: 
			if node is self: 
				continue
			dist_to_node = self.get_dist(node)
			no_NN = (self.dist_to_NN is None)
			if no_NN or (dist_to_node < self.dist_to_NN): 
				self.dist_to_NN = dist_to_node
				self.nn = node 
		return self.nn

# implementing UPGMA 
class UPGMA(object): 
	def __init__ (self, seqs, dist): 
		clusters = set([Node(dist, seq) for seq in seqs])
		self.dist = dist
		self.build_tree(clusters)
	
	#def getlargest(self): 
	#	if self.tree.right is None: 
	#		return self.tree.left 
	#	if self.tree.left is None: 
	#		return self.tree.right 
	#	elif len(self.tree.right)>len(self.tree.left): 
	#		return self.tree.right
	#	else: 
	#		return self.tree.left

	# Build UPGMA tree 
	def build_tree(self, clusters): 
		for node in clusters: 
			node.update_distances(clusters)

		# pick the closest pair and unionize them together 
		for i in range(len(clusters)-1): 
			c1,c2 = self.getClosestPair(clusters)
			#print("Clusters {} and {} closest".format(c1,c2))
			clusters.remove(c1)
			clusters.remove(c2)
			new_node = self.create_parentNode(c1,c2)
			clusters.add(new_node)
			new_node.update_distances(clusters)

			for node in clusters: 
				if (node.nn ==c1) or (node.nn==c2): 
					node.update_distances(clusters)
		assert len(clusters)==1 
		self.tree=clusters.pop()



	def create_parentNode (self,c1,c2): 
		return Node(self.dist, c1,c2)

	def getClosestPair(self, clusters): 
		min = None
		for node in clusters:
			if node.dist_to_NN == None: 
				node.dist_to_NN = 0

			if(min is None) or (node.dist_to_NN < min): 
				min = node.dist_to_NN
				c1=node
				c2=node.nn
		#print("Clusters {} and {} closest".format(c1,c2))
				#assert c2 is not None 
		return (c1, c2) 

	def findScore(c1,c2):
		match,mismatch, gap,extend = 0, 0, 0, 0
		gapMemory, subScore, counts1, counts2 = 0, 0, 0, 0
		#len1 = len(s1) - alignment.newLength(s1) #removes 'X'
		#len2 = len(s2) - alignment.newLength(s2) #removes 'X'
		#print(type(s2))
		#print(s2)
		len1 = len(s1)
		len2 = len(s2)
		for i in range(len1):
		    if (s1[i] or s2[i] == '-'):
		        if(s1[i] != '-' and s2[i] != '-'):
		            if (s1[i] == s2[i]): 
		                match += 1
		            else: 
		                mismatch += 1
		            counts1 = counts1 +1
		            counts2 = counts2 +1 
		            gapMemory=0      
		        elif (s1[i] != '-' and s2[i] == '-'):
		            if (counts2 >= 0 and counts2 < len2):
		                if(gapMemory == 1):
		                    extend = extend+1
		                else: 
		                    gap=gap+1
		                    gapMemory=1
		                counts1= counts1+1
		        elif (s1[i] == '-' and s2[i] != '-'):
		            if(counts1 >= 0 and counts1 < len1):
		                if(gapMemory == 2):
		                    extend = extend+1
		                else: 
		                    gap = gap +1
		                    gapMemory =2 
		                counts2 = counts2+1                    
		    subScore = (match*1) + (mismatch*(-1)) + (gap*(-4)) + (extend*(-1))
		    return subScore
	
	def totalScore(c1,c2):
		score = 0
		for node in c1: 
			for node2 in c2: 
				score += findScore(node, node2)
		return score 




if __name__=='__main__': 
	seqs = ['ACGT', 'GGTC', 'GTCA', 'GGCT']

	def scoringmatrix(s1,s2): 
		if (s1==s2): 
			return +1 
		if (s1!=s2): 
			return -1

	# gets distance between pairs
	def difference(node1, node2): 
		m = len(node1)
		n = len(node2)
		score = 0
		for x in range(m): 
			score += scoringmatrix(node1[x], node2[x])
		return score

	print(len(seqs))
	x=(UPGMA(seqs, difference).tree)
	# scoring for fitch gives you the best score between all the nodes 
	print(x)
	
	






