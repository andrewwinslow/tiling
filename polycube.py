
import unittest
import polyomino
import isohedral
import copy
import numpy
import sys
import math

# Input: an undirected graph represented as an adjacency set in a dict of sets,
#        and a vertex of the graph v.
# Output: whether the graph has any cycles in the component containing v.
def has_cycle(G, v):
	if len(G) == 0:
		return False
	visited = set([v])
	path = [v]
	def recurse(parent):
		for neigh in G[path[-1]]:
			if neigh == parent:
				continue
			if neigh in visited:
				return True
			path.append(neigh)
			visited.add(neigh)
			if recurse(path[-2]):
				return True
			path.pop()
		return False
	if recurse(None):
		return True
	return False
	

# Input: a directed graph represented as an adjacency set in a dict of sets
# Output: whether the graph is strongly connected
def is_connected(G):
	if len(G) == 0:
		return True
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

# Input: a set of edges represented as a set of 2-tuples.
# Output: an adjacency set representation of the undirected graph they induce.
def edge_set2graph(E):
	G = {}
	for v in [e[0] for e in E] + [e[1] for e in E]:
		G[v] = set([])
	for e in E:
		G[e[0]].add(e[1])
		G[e[1]].add(e[0])
	return G

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
# Output: The next clockwise half-edge around each face of the polycube's surface
def faces(P):
	incr_vecs = {'X': [(0, 1, 0), (0, 0, 1), (0, 1, 1)],
		'Y': [(1, 0, 0), (0, 0, 1), (1, 0, 1)],
		'Z': [(1, 0, 0), (0, 1, 0), (1, 1, 0)]}

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

def face_dual(P):

	def face_edges(f):
		bad_coords = filter(lambda c: math.floor(f[c] + 0.6) == f[c], [0, 1, 2])
		assert len(bad_coords) == 2
		okp = [list(f), list(f), list(f), list(f)]
		delta = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
		for i in xrange(4):
			okp[i][bad_coords[0]] = okp[i][bad_coords[0]] + 0.5*delta[i][0]
			okp[i][bad_coords[1]] = okp[i][bad_coords[1]] + 0.5*delta[i][1]
		return set([(tuple(okp[i]), tuple(okp[i+1])) for i in xrange(-1, 3)] + 
			[(tuple(okp[i+1]), tuple(okp[i])) for i in xrange(-1, 3)])

	V = set([])
	E = set([])
	unit_vecs = set([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
	for c in P:
		for vec in unit_vecs:
			for sign in [-1, 1]:
				adj_vec = tuple([sign*vec[i] for i in xrange(3)])
				if tuple([c[i] + adj_vec[i] for i in xrange(3)]) in P:
					continue
				vert = tuple([c[i] + 0.5*adj_vec[i] for i in xrange(3)])
				for v in V:
					if len(face_edges(vert) & face_edges(v)) > 0:
						E.add((vert, v))
				V.add(vert) 
	return E

	
# Input: a polycube represented as a set of integer 3-tuples
# Output: a generator of the boundary words of the polycube's unfoldings 
def unfoldings(P):
	faces_C = faces(P)
	faces_E = face_dual(P)
	faces_V = set([e[0] for e in faces_E] + [e[1] for e in faces_E])

	# Constructs a boundary word from a skeleton graph, 
	# randomly picking a component to trace the boundary of.	
	def skele_cycle_to_boundary_word(skele_cycle):
		if len(skele_cycle) == 0:
			return []

		def next_edge(cur_e, cur_d):	
			cur_d = polyomino.cw[cur_d]	
			naf = faces_C[cur_e]
			while not naf in skele_cycle:
				naf = faces_C[(naf[1], naf[0])] 
				cur_d = polyomino.ccw[cur_d]
			cur_e = naf	
			return cur_e, cur_d

		# Pick a starting edge and do an in-order traversal to define an Eulerian tour.
		start_e = list(skele_cycle)[0]
		W = ['N']
		cur_e, cur_dir = next_edge(start_e, 'N') 
		while not cur_e == start_e:
			W.append(cur_dir)	
			cur_e, cur_dir = next_edge(cur_e, cur_dir)
		return W

	def face_edges(f):
		# Construct all edges with this midpoint
		bad_coords = filter(lambda c: math.floor(f[c] + 0.6) == f[c], [0, 1, 2])
		okp = [list(f), list(f), list(f), list(f)]
		delta = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
		for i in xrange(4):
			okp[i][bad_coords[0]] = okp[i][bad_coords[0]] + 0.5*delta[i][0]
			okp[i][bad_coords[1]] = okp[i][bad_coords[1]] + 0.5*delta[i][1]
		return set([(tuple(okp[i]), tuple(okp[i+1])) for i in xrange(-1, 3)] + 
			[(tuple(okp[i+1]), tuple(okp[i])) for i in xrange(-1, 3)])

	def graph_to_boundary_word(T):
		T_V = set([e[0] for e in T] + [e[1] for e in T])
		# Take dual edges of those not in T
		comp_T = faces_E - T 
		# Construct the cycle through these edges that traces
		# the boundary of T clockwise
		skele_cycle = set([])
		for e in comp_T:
			# If edge doesn't connect to something in T
			if (not e[0] in T_V) and (not e[1] in T_V):
				continue 
			cand_edges = list(set(face_edges(e[0]) & face_edges(e[1]) & set(faces_C.keys())))	
			if faces_C[cand_edges[0]] in face_edges(e[1]):
				cand_edges = [cand_edges[1], cand_edges[0]] 
			if e[0] in T_V:
				skele_cycle.add(cand_edges[0])
			if e[1] in T_V:
				skele_cycle.add(cand_edges[1])
		return skele_cycle_to_boundary_word(skele_cycle)	

	pT = set([])	
	rem_E = sorted(list(faces_E))
	G_pT = {}
	for v in faces_V:
		G_pT[v] = set([])
	G_T_rem_E = {}
	for v in faces_V:
		G_T_rem_E[v] = set([])
		for e in faces_E:
			if v == e[0]:
				G_T_rem_E[v].add(e[1])
			if v == e[1]:
				G_T_rem_E[v].add(e[0])

	def nonplanar_pT():
		return not polyomino.is_weakly_simple(graph_to_boundary_word(pT))
		
	def recurse():
		if len(pT) == len(faces_V) - 1:
			yield pT
			return
		if len(rem_E) == 0:
			return
		# Adversarily decide on an edge b to branch on.
		# Thus look for one that kills one of the two branches.
		b = '?'
		for pb in rem_E:
			G_T_rem_E[pb[0]].remove(pb[1])
			G_T_rem_E[pb[1]].remove(pb[0])
			if not is_connected(G_T_rem_E):
				b = pb
			G_T_rem_E[pb[0]].add(pb[1])
			G_T_rem_E[pb[1]].add(pb[0])
			
			G_pT[pb[0]].add(pb[1])
			G_pT[pb[1]].add(pb[0])
			if (has_cycle(G_pT, pb[0]) or nonplanar_pT()):
				b = pb	
			G_pT[pb[0]].remove(pb[1])
			G_pT[pb[1]].remove(pb[0])
			if b != '?':
				break
		if b == '?':
			b = rem_E[0]
		rem_E.remove(b)

		# Recursion 1: b is not included. 
		# Recurse with slightly smaller remaining edge set.
		G_T_rem_E[b[0]].remove(b[1])
		G_T_rem_E[b[1]].remove(b[0])
		if is_connected(G_T_rem_E):
			for T in recurse():
				yield T
		G_T_rem_E[b[0]].add(b[1])
		G_T_rem_E[b[1]].add(b[0])
		# Recursion 2: b is included.
		# Recurse with slightly smaller remaining edge set and tree
		pT.add(b)
		G_pT[b[0]].add(b[1])
		G_pT[b[1]].add(b[0])
		# Gut check: does T have any cycles? Die if so.
		# For speed, only check component containing b
		if not (has_cycle(G_pT, b[0]) or nonplanar_pT()): 
			for T in recurse():
				yield T
		# Restore variables
		pT.remove(b)
		G_pT[b[0]].remove(b[1])
		G_pT[b[1]].remove(b[0])
		rem_E.append(b)

	for T in recurse():
		yield graph_to_boundary_word(T)


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
		count = 0
		for uf in unfoldings(set([(0, 0, 0)])):
			count = count + 1
			self.assertTrue(polyomino.is_polyomino(uf))
			self.assertTrue(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))
		self.assertEqual(count, ((2*1)**3 * (2*2)**3 * (2*3)**1) / 8) 
		# Some unfolding of a tricube
		tiler = False
		for uf in unfoldings(set([(0, 0, 0), (1, 0, 0), (2, 0, 0)])):
			cuf = polyomino.cancel(uf)
			if (polyomino.is_polyomino(cuf) and 
				(isohedral.has_translation_tiling(cuf) or isohedral.has_half_turn_tiling(cuf))):
				tiler = True
				break
		self.assertTrue(tiler)

	def test__unfolding_count(self):
		self.assertEqual(unfolding_count(set([(0, 0, 0)])), ((2*1)**3 * (2*2)**3 * (2*3)**1) / 8)

	def test__face_dual(self):
		E = face_dual(set([(0, 0, 0)]))
		self.assertEqual(len(E), 12) 

	def test__has_cycle(self):
		G = {1: set([2]), 2: set([1, 3]), 3: set([2])}
		self.assertFalse(has_cycle(G, 1))
		self.assertFalse(has_cycle(G, 2))
		self.assertFalse(has_cycle(G, 3))
		G[3].add(1)
		G[1].add(3)
		self.assertTrue(has_cycle(G, 1))
		self.assertTrue(has_cycle(G, 2))
		self.assertTrue(has_cycle(G, 3))

		G = {1: set([2]), 2: set([1, 3]), 3: set([2, 4]), 4: set([3])}
		self.assertFalse(has_cycle(G, 1))
		self.assertFalse(has_cycle(G, 2))
		self.assertFalse(has_cycle(G, 3))
		G[4].add(1)
		G[1].add(4)
		self.assertTrue(has_cycle(G, 1))
		self.assertTrue(has_cycle(G, 2))
		self.assertTrue(has_cycle(G, 3))
		self.assertTrue(has_cycle(G, 4))

		self.assertFalse(has_cycle({1: set([2]), 2: set([1])}, 1))
		self.assertFalse(has_cycle({1: set([2]), 2: set([1])}, 2))

if __name__ == '__main__':
	unittest.main()
	#print len([uf for uf in unfoldings(set([(0, 0, 0), (1, 0, 0)]))])

