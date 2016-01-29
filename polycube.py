
import unittest
import polyomino
import isohedral
import copy
import numpy

# Input: a directed graph represented as an adjacency list in a dict of sets
# Output: whether the graph is strongly connected
def is_connected(G):
	root = G.keys()[0]
	visited = set([root])
	path = [root]
	def recurse():
		for neigh in G[path[-1]]:
			if not (neigh in visited):
				path.append(neigh)
				visited.add(neigh)
				recurse()
				path.pop()
	recurse()
	return len(visited) == len(G)		

# Input: a set of integer 3-tuples
# Output: whether the set of cells their describe has a connected dual graph
def is_polycube(P):
	def neighbors(cell):
		vecs = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
		return set([(cell[0] + v[0], cell[1] + v[1], cell[2] + v[2]) for v in vecs])

	# Construct the dual graph
	G = {}
	for cell in P:
		G[cell] = set([])
	for cell in P:
		for neigh in neighbors(cell) & P:
			G[cell].add(neigh)
			G[neigh].add(cell)
	return is_connected(G)

# Input: a polycube represented as a set of integer 3-tuples
# Output: The set of edges of the skeleton
def skeleton(P):
	incr_vecs = {'X': [(0, 1, 0), (0, 0, 1), (0, 1, 1)],
		'Y': [(1, 0, 0), (0, 0, 1), (1, 0, 1)],
		'Z': [(1, 0, 0), (0, 1, 0), (1, 1, 0)]}
	def neighbors(v, d):
		return set([(v[0]+vec[0], v[1]+vec[1], v[2]+vec[2]) for vec in incr_vecs[d] + [(0, 0, 0)]])

	E = set([])
	x_range = (min([x for (x, y, z) in P]), max([x for (x, y, z) in P]))
	y_range = (min([y for (x, y, z) in P]), max([y for (x, y, z) in P]))
	z_range = (min([z for (x, y, z) in P]), max([z for (x, y, z) in P]))
	for x in xrange(x_range[0]-1, x_range[1]+1):
		for y in xrange(y_range[0]-1, y_range[1]+1):
			for z in xrange(z_range[0]-1, z_range[1]+1):
				# For each edge, check if it's incident to between 1 and 3 cells
				# If so, add it to the skeleton.
				if 0 < len(neighbors((x, y, z), 'X') & P) < 4:
					E.add(((x-0.5, y+0.5, z+0.5), (x+0.5, y+0.5, z+0.5)))	
				if 0 < len(neighbors((x, y, z), 'Y') & P) < 4:
					E.add(((x+0.5, y-0.5, z+0.5), (x+0.5, y+0.5, z+0.5)))	
				if 0 < len(neighbors((x, y, z), 'Z') & P) < 4:
					E.add(((x+0.5, y+0.5, z-0.5), (x+0.5, y+0.5, z+0.5)))	
	return E

# Input: a polycube represented as a set of integer 3-tuples
# Output: the number of unfoldings
def unfolding_count(P):
	def determinant(A):
		return int(numpy.linalg.det(A))
		"""
		# Just takes wayyyyy too long...
		if len(A) == 2:
			return A[0][0]*A[1][1] - A[0][1]*A[1][0]
		total = 0
		for i in xrange(len(A)):
			factor = A[0][i]*(-1)**(i % 2)
			if factor != 0: # A speedup filter
				new_A = copy.deepcopy(A)
				del new_A[0]
				for j in xrange(len(new_A)):
					del new_A[j][i]
				total = total + factor*determinant(new_A)
		return total
		"""

	skele_E = skeleton(P)	
	skele_V = list(set([e[0] for e in skele_E] + [e[1] for e in skele_E]))
	A = [[0] * len(skele_V) for i in xrange(len(skele_V))]
	for i in xrange(len(skele_V)):
		for j in xrange(len(skele_V)):
			e = (skele_V[i], skele_V[j])
			if e in skele_E or (e[1], e[0]) in skele_E:
				A[i][j] = -1 # Invert now, since we really want -A
	# Set diagonals to degrees (negative sum of -1 entries in row)
	for i in xrange(len(skele_V)):
		A[i][i] = -sum(A[i])	
	# Now delete one row and column (the last ones)
	for i in xrange(len(skele_V)):
		del A[i][-1]
	del A[-1]	
	# Now return determinant
	return determinant(A)
	

# Input: a polycube represented as a set of integer 3-tuples
# Output: The next half-edge around each face of the polycube's surface
def faces(P):
	incr_vecs = {'X': [(0, 1, 0), (0, 0, 1), (0, 1, 1)],
		'Y': [(1, 0, 0), (0, 0, 1), (1, 0, 1)],
		'Z': [(1, 0, 0), (0, 1, 0), (1, 1, 0)]}
	def neighbors(v, d):
		return set([(v[0]+vec[0], v[1]+vec[1], v[2]+vec[2]) for vec in incr_vecs[d] + [(0, 0, 0)]])

	G = {}
	for c in P:
		if not (c[0] + 1, c[1], c[2]) in P:
			verts = [(c[0]+0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]-0.5, c[2]-0.5),
				(c[0]+0.5, c[1]-0.5, c[2]+0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0] - 1, c[1], c[2]) in P:
			verts = [(c[0]-0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]-0.5, c[2]+0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]+0.5, c[2]-0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1]+1, c[2]) in P:
			verts = [(c[0]-0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]+0.5, c[2]+0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1]-1, c[2]) in P:
			verts = [(c[0]+0.5, c[1]-0.5, c[2]+0.5), 
				(c[0]+0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]+0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1], c[2]+1) in P:
			verts = [(c[0]-0.5, c[1]-0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]+0.5, c[2]+0.5),
				(c[0]+0.5, c[1]+0.5, c[2]+0.5),
				(c[0]+0.5, c[1]-0.5, c[2]+0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1], c[2]-1) in P:
			verts = [(c[0]+0.5, c[1]-0.5, c[2]-0.5), 
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]-0.5, c[1]+0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5)]
			for i in xrange(4):
				G[(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
	return G






# Input: a polycube represented as a set of integer 3-tuples
# Output: a generator of the boundary words of the polycube's unfoldings 
def unfoldings(P):
	skele_E = skeleton(P)
	faces_G = faces(P)
	skele_V = set([e[0] for e in skele_E] + [e[1] for e in skele_E])
	n = len(skele_V)

	# Label the edges of the skeleton
	edge2id = {}
	id2edge = {}
	for e in skele_E:
		eid = len(edge2id)
		edge2id[e] = eid
		id2edge[eid] = e

	# Enumerate the spanning trees using the Minty algorithm
	# from G. J. Minty, "A simple algorithm for listing all the trees of a graph",
	# IEEE Transactions on Circuit Theory, 12(1), pp. 120, 1965.
	edge_ids = []
	def enumerate_trees(E, edge2id):
		if len(E) == 0:
			yield set([id2edge[eid] for eid in edge_ids])
			return
		# Choose an edge b. Will recurse in two ways:
		# 1. Build a graph G1 by contracting b (adding it to spanning tree)
		# 2. Build a graph G2 by deleting b (not adding it to spanning tree)
		b = list(E)[0]
		bid = len(edge_ids)
		# G1: shrink e to a vertex, add e to partial edge set, recurse
		new_E = set([])
		new_edge2id = {}
		rlbl = {}
		for v in set([e[0] for e in E] + [e[1] for e in E]):
			rlbl[v] = v
		rlbl[b[0]] = bid
		rlbl[b[1]] = bid
		for e in E:
			new_e = (rlbl[e[0]], rlbl[e[1]])
			if new_e[0] != new_e[1]:
				new_E.add(new_e)
				new_edge2id[new_e] = edge2id[e]
		edge_ids.append(edge2id[b])
		for T in enumerate_trees(new_E, new_edge2id):
			yield T
		edge_ids.pop()
		# G2: delete e, check if remainder is connected, recurse if so
		new_E = E.copy()
		new_E.remove(b)
		new_edge2id = edge2id.copy()
		del new_edge2id[b]
		new_G = {}
		for V in set([e[0] for e in E] + [e[1] for e in E]):
			new_G[V] = set([])
		for e in new_E:
			new_G[e[0]].add(e[1])
			new_G[e[1]].add(e[0])
		if is_connected(new_G):
			for T in enumerate_trees(new_E, new_edge2id):
				yield T

	# Iterate over subtrees
	for tree in enumerate_trees(skele_E, edge2id):
		half_edges = set([(e[0], e[1]) for e in tree] + [(e[1], e[0]) for e in tree])
		rot = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}
		def next_edge(cur_e, cur_d):	
			while not (faces_G[cur_e] in half_edges):
				along_face = faces_G[cur_e] 
				cur_e = (along_face[1], along_face[0])
				cur_d = rot[rot[rot[cur_d]]]
			cur_e = faces_G[cur_e]
			cur_d = rot[cur_d]
			return cur_e, cur_d

		# Pick a starting edge and do an in-order traversal to define an Eulerian tour
		cur_e = list(half_edges)[0]
		cur_dir = 'N'
		W = []
		for i in xrange(2*(n-1)):
			cur_e, cur_dir = next_edge(cur_e, cur_dir)
			W.append(cur_dir)	
		polyomino.cancel(W)
		yield W


class TestStuff(unittest.TestCase):
	
	def setUp(self):
		pass

	def test__is_polycube(self):
		self.assertTrue(is_polycube(set([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])))
		self.assertTrue(is_polycube(set([(1, 1, 1), (1, 1, 2), (1, 2, 2)])))
		self.assertFalse(is_polycube(set([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)])))

	def test__skeleton(self):	
		S = skeleton(set([(0, 0, 0)])) 
		self.assertEqual(len(S), 12)
		h = 0.5
		p = [(-h, -h, -h), (-h, h, -h), (h, h, -h), (h, -h, -h), (-h, -h, h), (-h, h, h), (h, h, h), (h, -h, h)]
		# Bottom edges
		self.assertTrue((p[0], p[1]) in S or (p[1], p[0]) in S)
		self.assertTrue((p[1], p[2]) in S or (p[2], p[1]) in S)
		self.assertTrue((p[2], p[3]) in S or (p[3], p[2]) in S)
		self.assertTrue((p[0], p[3]) in S or (p[3], p[0]) in S)
		# Top edges
		self.assertTrue((p[4], p[5]) in S or (p[5], p[4]) in S)
		self.assertTrue((p[5], p[6]) in S or (p[6], p[5]) in S)
		self.assertTrue((p[6], p[7]) in S or (p[7], p[6]) in S)
		self.assertTrue((p[4], p[7]) in S or (p[7], p[4]) in S)
		# Vertical edges
		self.assertTrue((p[0], p[4]) in S or (p[4], p[0]) in S)
		self.assertTrue((p[1], p[5]) in S or (p[5], p[1]) in S)
		self.assertTrue((p[2], p[6]) in S or (p[6], p[2]) in S)
		self.assertTrue((p[3], p[7]) in S or (p[7], p[3]) in S)
			
		S = skeleton(set([(0, 0, 0), (0, 0, 1)]))
		th = 1.5
		p = [(-h, -h, -h), (-h, h, -h), (h, h, -h), (h, -h, -h), 
			(-h, -h, h), (-h, h, h), (h, h, h), (h, -h, h),
			(-h, -h, th), (-h, h, th), (h, h, th), (h, -h, th)]
		# Bottom edges
		self.assertTrue((p[0], p[1]) in S or (p[1], p[0]) in S)
		self.assertTrue((p[1], p[2]) in S or (p[2], p[1]) in S)
		self.assertTrue((p[2], p[3]) in S or (p[3], p[2]) in S)
		self.assertTrue((p[0], p[3]) in S or (p[3], p[0]) in S)
		# Middle edges
		self.assertTrue((p[4], p[5]) in S or (p[5], p[4]) in S)
		self.assertTrue((p[5], p[6]) in S or (p[6], p[5]) in S)
		self.assertTrue((p[6], p[7]) in S or (p[7], p[6]) in S)
		self.assertTrue((p[4], p[7]) in S or (p[7], p[4]) in S)
		# Top edges
		self.assertTrue((p[8], p[9]) in S or (p[8], p[9]) in S)
		self.assertTrue((p[9], p[10]) in S or (p[10], p[9]) in S)
		self.assertTrue((p[10], p[11]) in S or (p[11], p[10]) in S)
		self.assertTrue((p[8], p[11]) in S or (p[11], p[8]) in S)
		# Lower vertical edges
		self.assertTrue((p[0], p[4]) in S or (p[4], p[0]) in S)
		self.assertTrue((p[1], p[5]) in S or (p[5], p[1]) in S)
		self.assertTrue((p[2], p[6]) in S or (p[6], p[2]) in S)
		self.assertTrue((p[3], p[7]) in S or (p[7], p[3]) in S)
		# Upper vertical edges
		self.assertTrue((p[4], p[8]) in S or (p[8], p[4]) in S)
		self.assertTrue((p[5], p[9]) in S or (p[9], p[5]) in S)
		self.assertTrue((p[6], p[10]) in S or (p[10], p[6]) in S)
		self.assertTrue((p[7], p[11]) in S or (p[11], p[7]) in S)

	def test__unfoldings(self):
		# All unfoldings of cube tile
		for uf in unfoldings(set([(0, 0, 0)])):
			self.assertTrue(polyomino.is_polyomino(uf))
			self.assertTrue(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))
		# Some unfolding of dicube tiles
		tiler = False
		for uf in unfoldings(set([(0, 0, 0), (0, 0, 1)])):
			if (polyomino.is_polyomino(uf) and 
				(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))):
				tiler = True
				break
		self.assertTrue(tiler)
		# Some unfolding of another dicube tiles
		tiler = False
		for uf in unfoldings(set([(0, 0, 0), (1, 0, 0)])):
			if (polyomino.is_polyomino(uf) and 
				(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))):
				tiler = True
				break
		self.assertTrue(tiler)
		# Some unfolding of a tricube
		tiler = False
		for uf in unfoldings(set([(0, 0, 0), (1, 0, 0), (0, 1, 0)])):
			if (polyomino.is_polyomino(uf) and 
				(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))):
				tiler = True
				break
		self.assertTrue(tiler)
		# Some unfolding of a tricube
		tiler = False
		for uf in unfoldings(set([(0, 0, 0), (1, 0, 0), (2, 0, 0)])):
			if (polyomino.is_polyomino(uf) and 
				(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))):
				tiler = True
				break
		self.assertTrue(tiler)

	def test__unfolding_count(self):
		self.assertEqual(unfolding_count(set([(0, 0, 0)])), ((2*1)**3 * (2*2)**3 * (2*3)**1) / 8)

if __name__ == '__main__':
	unittest.main()

