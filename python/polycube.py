
import unittest
import polyomino
import isohedral
import copy
try:
	import numpy
except Exception:
	numpy = None
import sys
import math
import time
import random

# Input: an undirected graph represented as an adjacency set in a dict of sets,
#        and a vertex of the graph v.
# Output: whether the graph has any cycles in the component containing v.
def has_cycle(G, v):
	if len(G) == 0:
		return False
	visited = []
	path = [v]
	def recurse(parent):
		for neigh in G[path[-1]]:
			if neigh == parent:
				continue
			if neigh in visited:
				return True
			path.append(neigh)
			visited.append(neigh)
			if recurse(path[-2]):
				return True
			path.pop()
		return False
	if recurse(None):
		return True
	return False
	
# Input: a directed graph represented as an adjacency set in a dict of sets
#        and a vertex of this graph. 
# Output: the number of vertices reachable from the vertex.
def reachable(G, root):
	visited = [root]
	path = [root]
	def recurse():
		for neigh in G[path[-1]]:
			if not (neigh in visited):
				path.append(neigh)
				visited.append(neigh)
				recurse()
				path.pop()
	recurse()
	return visited

# Input: a directed graph represented as an adjacency set in a dict of sets
# Output: whether the graph is strongly connected
def is_connected(G):
	if len(G) == 0:
		return True
	root = G.keys()[0]
	visited = [root]
	path = [root]
	def recurse():
		for neigh in G[path[-1]]:
			if not (neigh in visited):
				path.append(neigh)
				visited.append(neigh)
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

# Input: a polycube described by a set of integer 3-tuples
# Output: the cell dual graph
def cell_dual(P):
	def neighbors(cell, P):
		neighs = []
		for vec in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
			adj = (cell[0] + vec[0], cell[1] + vec[1], cell[2] + vec[2]) 
			if adj in P:
				neighs.append(adj)
		return neighs

	# Construct the dual graph
	G = {}
	for cell in P:
		G[cell] = set([])
	for cell in P:
		for neigh in neighbors(cell, P):
			G[cell].add(neigh)
			G[neigh].add(cell)
	return G

# Input: a set of integer 3-tuples
# Output: whether the set of cells their describe has a connected dual graph
def is_polycube(P):
	return is_connected(cell_dual(P))

# Input: an undirected graph represented as an dict of adjacency sets
# Output: the number of spanning trees of the graph
def spanning_tree_count(G):
	
	def determinant(A):
		if numpy:
			return int(numpy.linalg.det(A))
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
		
	V = list(G.keys()) # Need V to be indexable
	if len(V) <= 1:
		return 1
	A = [[0] * len(V) for i in xrange(len(V))]
	for i in xrange(len(V)):
		for j in xrange(len(V)):
			if V[j] in G[V[i]]:
				A[i][j] = -1 # Invert now, since we really want -A
	# Set diagonals to degrees (negative sum of -1 entries in row)
	for i in xrange(len(V)):
		A[i][i] = -sum(A[i])	
	# Now delete one row and column (the last ones)
	for i in xrange(len(V)):
		del A[i][-1]
	del A[-1]	
	# Now return determinant
	return determinant(A)

# Input: a polycube represented as a set of integer 3-tuples
# Output: the number of unfoldings
def unfolding_count(P):
	V, E = face_dual(P)
	G = {}
	for v in V:
		G[v] = set([])
		for e in E:
			if e[0] == v:
				G[v].add(e[1])
			elif e[1] == v:
				G[v].add(e[0])
	return spanning_tree_count(G)	

# Input: a polycube represented as a set of integer 3-tuples
# Output: Two dictionaries. The first has a set of clockwise half-edges for each face and
#         the next clockwise edge for each edge. The second dictionary is the same, 
#         but for counterclockwise edges.
def face_cycles(P):
	G = {}
	for c in P:
		if not (c[0] + 1, c[1], c[2]) in P:
			verts = [(c[0]+0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]-0.5, c[2]-0.5),
				(c[0]+0.5, c[1]-0.5, c[2]+0.5)]
			f = (c[0]+0.5, c[1], c[2])
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])

		if not (c[0] - 1, c[1], c[2]) in P:
			verts = [(c[0]-0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]-0.5, c[2]+0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]+0.5, c[2]-0.5)]
			f = (c[0]-0.5, c[1], c[2])
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1]+1, c[2]) in P:
			verts = [(c[0]-0.5, c[1]+0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]+0.5, c[1]+0.5, c[2]+0.5)]
			f = (c[0], c[1]+0.5, c[2])
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1]-1, c[2]) in P:
			verts = [(c[0]+0.5, c[1]-0.5, c[2]+0.5), 
				(c[0]+0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]+0.5)]
			f = (c[0], c[1]-0.5, c[2])
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1], c[2]+1) in P:
			verts = [(c[0]-0.5, c[1]-0.5, c[2]+0.5), 
				(c[0]-0.5, c[1]+0.5, c[2]+0.5),
				(c[0]+0.5, c[1]+0.5, c[2]+0.5),
				(c[0]+0.5, c[1]-0.5, c[2]+0.5)]
			f = (c[0], c[1], c[2]+0.5)
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
		if not (c[0], c[1], c[2]-1) in P:
			verts = [(c[0]+0.5, c[1]-0.5, c[2]-0.5), 
				(c[0]+0.5, c[1]+0.5, c[2]-0.5),
				(c[0]-0.5, c[1]+0.5, c[2]-0.5),
				(c[0]-0.5, c[1]-0.5, c[2]-0.5)]
			f = (c[0], c[1], c[2]-0.5)
			G[f] = {}
			for i in xrange(4):
				G[f][(verts[i-2], verts[i-1])] = (verts[i-1], verts[i])
	CCW_G = {}
	for f in G:
		CCW_G[f] = {}
		for e in G[f]:
			next_e = G[f][e]
			CCW_G[f][(next_e[1], next_e[0])] = (e[1], e[0])
	return G, CCW_G

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

	cell_G = cell_dual(P)
	def cell_distance_leq_2(c1, c2):
		if c1 == c2:
			return True
		if c2 in cell_G[c1]:
			return True
		for c3 in cell_G[c1]:
			if c2 in cell_G[c3]:
				return True
		return False

	v2c = {}
	V = set([])
	E = set([])
	for c in P:
		for vec in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
			for sign in [-1, 1]:
				adj_vec = tuple([sign*vec[i] for i in [0, 1, 2]])
				if tuple([c[i] + adj_vec[i] for i in [0, 1, 2]]) in P:
					continue
				vert = tuple([c[i] + 0.5*adj_vec[i] for i in [0, 1, 2]])
				v2c[vert] = c
				for v in V:
					# Must share a common edge and not do so degenerately
					# Non-degenerate: the faces' cells are either the same (reflex), 
					# adjacent (flat), or distance 2 (convex)  
					if (len(face_edges(vert) & face_edges(v)) > 0 and 
						cell_distance_leq_2(v2c[vert], v2c[v])):
						E.add((vert, v))
				V.add(vert) 
	return V, E

# Input: a polycube and a weakly simple boundary word
# Output: whether the weakly simple polyomino given is an unfolding
#         of the polycube's surface.
def is_unfolding(P, U):
	faces_CW, faces_CCW = face_cycles(P)
	P_dual_V, P_dual_E = face_dual(P)

	P_G = {}
	for v in P_dual_V:
		P_G[v] = []
	for e in P_dual_E:
		P_G[e[0]].append(e[1])
		P_G[e[1]].append(e[0])
	
	# Fix a root of the face dual
	max_z = max([p[2] for p in P_G.keys()])
	P_root = min(filter(lambda p: p[2] == max_z, P_G.keys()))	

	def consistent(P_G, P_root, U_G, U_root):
		# Do a syncronized DFS through both graphs
		P_visited = [P_root]
		U_visited = [U_root]
		P_path = [P_root]
		U_path = [U_root]
		def recurse(cur_P_d, cur_P_e):
			for U_neigh in U_G[U_path[-1]]:
				if U_neigh in U_visited:
					continue
				next_P_d = cur_P_d
				next_P_e = cur_P_e
				# Determine direction of neigh
				neigh_vec = (U_neigh[0] - U_path[-1][0], U_neigh[1] - U_path[-1][1]) 
				neigh_d = polyomino.vec2dir[neigh_vec]
				# Compute adjacent face in P of same dir
				while next_P_d != neigh_d:
					next_P_e = faces_CW[P_path[-1]][next_P_e]
					next_P_d = polyomino.cw[next_P_d]
				P_neigh = '?'
				for v in P_G[P_path[-1]]:
					if next_P_e in faces_CCW[v]:
						P_neigh = v
						next_P_d = polyomino.comp[next_P_d]
						next_P_e = (next_P_e[1], next_P_e[0])
						break
				assert P_neigh != '?'
				U_path.append(U_neigh)
				U_visited.append(U_neigh)
				P_path.append(P_neigh)
				P_visited.append(P_neigh)
				recurse(next_P_d, next_P_e)
				U_path.pop()
				P_path.pop()
		# Start with arbitrary choice of north meaning a specific edge of the root
		# Note: must pick this edge the same way for all inputs with same P_root
		recurse('N', min(faces_CW[P_root]))
		return len(set(P_visited)) == len(P_G) 

	# Assumes north is the same in P and U
	# Returns the cell that's mapped to the face root or False
	def has_unfolding(P, P_root, U):
		U_dual_V, U_dual_E = polyomino.cell_dual(U)
		U_G = {}
		for v in U_dual_V:
			U_G[v] = []
		for e in U_dual_E:
			U_G[e[0]].append(e[1])	
			U_G[e[1]].append(e[0])	

		for U_root in U_G:
			if consistent(P_G, P_root, U_G, U_root):
				return U_root
		return False
		
	# Iterate through all possible orientations
	for rot in [0, 90, 180, 270]:
		rotU = [polyomino.rot[rot][d] for d in U]
		result = has_unfolding(P, P_root, rotU)
		if result:
			return (P_root, result, polyomino.cell_dual(rotU)[0])
		reflrotU = [polyomino.refl[0][d] for d in rotU]
		result = has_unfolding(P, P_root, reflrotU) 
		if result:
			return (P_root, result, polyomino.cell_dual(reflrotU)[0])
		
	return False
	
# Input: a polycube represented as a set of integer 3-tuples
# Output: a generator of the boundary words of the polycube's unfoldings 
def unfoldings(P, strongly_simple=False, hamiltonian=False):
	faces_CW, faces_CCW = face_cycles(P)
	faces_V, faces_E = face_dual(P)

	rem_E = sorted(list(faces_E)) # Sort to bias towards "filling out" along X-axis?
	# Initialize graph of partial unfolding tree
	G_pT = {}
	for v in faces_V:
		G_pT[v] = set([])
	# Initialize face dual graph of rest of surface 
	G_T_rem_E = {}
	for v in faces_V:
		G_T_rem_E[v] = set([])
		for e in faces_E:
			if v == e[0]:
				G_T_rem_E[v].add(e[1])
			if v == e[1]:
				G_T_rem_E[v].add(e[0])
	
	def tree_to_boundary_word():
		if sum([len(v) for v in G_pT.values()]) < 2:
			return ['N', 'E', 'S', 'W']
	
		start_e = '?'
		cur_v = '?' 
		for v in G_pT:
			# Find a face edge of a vertex in pT that 
			# borders another face not adjacent in pT
			if len(G_pT[v]) != 1:
				continue
			start_v = v
			for e in faces_CW[v]:
				if not (e in faces_CCW[min(G_pT[v])]) and not (faces_CW[v][e] in faces_CCW[min(G_pT[v])]):
					start_e = e
			break	

		W = [polyomino.cw['S']]
		W_edges = [start_e]
		cur_v = start_v
		cur_d = polyomino.cw['S']
		cur_e = faces_CW[cur_v][start_e]
		while (cur_e, cur_v) != (start_e, start_v):
			edge2neigh_v = {}
			for neigh_v in G_pT[cur_v]:
				edge2neigh_v[min(set(faces_CCW[neigh_v]) & set(faces_CW[cur_v]))] = neigh_v	
			if cur_e in edge2neigh_v: # Not an edge of tree's boundary
				# Move to the adjoining cell
				cur_v = edge2neigh_v[cur_e]
				cur_e = faces_CW[cur_v][(cur_e[1], cur_e[0])]
				cur_d = polyomino.cw[polyomino.comp[cur_d]]
			else:	# An edge of the tree's boundary
				# Add the current edge to the boundary word
				W.append(polyomino.cw[cur_d])
				W_edges.append(cur_e)
				cur_d = polyomino.cw[cur_d]
				cur_e = faces_CW[cur_v][cur_e]
		# Glue coincident edges back together	
		welded = False
		while not welded:
			welded = True
			for i in xrange(len(W)-1):
				if W[i] == polyomino.comp[W[i+1]] and W_edges[i][0] == W_edges[i+1][1] and W_edges[i][1] == W_edges[i+1][0]:
					del W[i+1] 
					del W_edges[i+1]
					del W[i]
					del W_edges[i]
					welded = False
					break
			if W[0] == polyomino.comp[W[-1]] and W_edges[0][0] == W_edges[-1][1] and W_edges[0][1] == W_edges[-1][0]:
				del W[-1]
				del W_edges[-1]
				del W[0]
				del W_edges[0]
				welded = False
		return W

	def skipped_count():
		contract_G = {'pT': set([])}
		vmap = {}
		for v in G_T_rem_E:
			if len(G_pT[v]) > 0: # If in partial (connected) tree
				vmap[v] = 'pT'
			else:
				vmap[v] = v
				contract_G[v] = set([])
		for v in G_T_rem_E:
			for neigh in G_T_rem_E[v]:
				contract_G[vmap[v]].add(vmap[neigh])
		return spanning_tree_count(contract_G)	

	stats = [0, 0]
	def recurse():
		#sys.stdout.write("\rProgress: " + str(stats[0]+stats[1]) + " (" + str(stats[0]) + ") / " + str(total_unfoldings))
		#sys.stdout.flush()
		pT_edge_count = sum([len(G_pT[v]) for v in G_pT]) / 2
		# If the number of edges is right
		if pT_edge_count == len(faces_V) - 1:
			stats[0] = stats[0] + 1
			yield tree_to_boundary_word()
			return
		# If there aren't enough edges left to make a tree
		if len(rem_E) + pT_edge_count < len(faces_V) - 1:
			return
		# Adversarily decide on an edge b to branch on.
		# Look for one that kills one or both(!) of the two branches.
		branch1_killer = '?' # Exclusion is impossible
		branch2_killer = '?' # Inclusion is impossible

		# Assumes pb has already been added to G_pT
		def kills_branch2(pb):
			if has_cycle(G_pT, pb[0]):
				return True
			BW = tree_to_boundary_word()
			if strongly_simple and not polyomino.is_simple(BW):
				return True
			if not polyomino.is_weakly_simple(tree_to_boundary_word()):
				return True
			if hamiltonian:
				# If the partial tree is not a path
				if len(filter(lambda v: len(G_pT[v]) > 2, G_pT.keys())) > 0:
					return True
				# If the partial unfolding isn't a path
				# (two faces not adjacent on the path are adjacent on surface)
				if len(BW) != 2 + 2 * len(filter(lambda v: len(G_pT[v]) > 0, G_pT.keys())): 
					return True
			return False			

		for pb in rem_E:
			# These must necessarily be connected to existing partial tree
			if len(G_pT[pb[0]]) + len(G_pT[pb[1]]) == 0:
				continue
			G_T_rem_E[pb[0]].remove(pb[1])
			G_T_rem_E[pb[1]].remove(pb[0])
			if not is_connected(G_T_rem_E):
				branch1_killer = pb
			G_T_rem_E[pb[0]].add(pb[1])
			G_T_rem_E[pb[1]].add(pb[0])
			
			G_pT[pb[0]].add(pb[1])
			G_pT[pb[1]].add(pb[0])
			if kills_branch2(pb):
				branch2_killer = pb	
			G_pT[pb[0]].remove(pb[1])
			G_pT[pb[1]].remove(pb[0])

			# If you can kill both branches, do it
			if branch1_killer == pb and branch2_killer == pb:
				break
		
		# Take branch 2 killer preferably (branch 1 is closer to end b/c edge starvation)
		b = '?'
		if branch2_killer != '?':
			b = branch2_killer
		elif branch1_killer != '?':
			b = branch1_killer
		# If neither branch can be killed, pick something
		# that's connected to the partial tree that's growing
		if b == '?':
			if pT_edge_count == 0: 
				b = random.choice(rem_E) 
			else:
				b = random.choice(filter(lambda e: len(G_pT[e[0]]) + len(G_pT[e[1]]) > 0, rem_E))
		rem_E.remove(b)
	
		# Recursion branch 1: b is not included. 
		# Recurse with slightly smaller remaining edge set.
		G_T_rem_E[b[0]].remove(b[1])
		G_T_rem_E[b[1]].remove(b[0])
		# Branch killer 1: cannot possibly finish a (connected) tree?
		if is_connected(G_T_rem_E): 
			for W in recurse():
				yield W
		G_T_rem_E[b[0]].add(b[1])
		G_T_rem_E[b[1]].add(b[0])
		# Recursion branch 2: b is included.
		# Recurse with slightly smaller remaining edge set and tree
		G_pT[b[0]].add(b[1])
		G_pT[b[1]].add(b[0])
		# Branch killer 2: T has a cycles, a nonplanar unfolding, 
		# or is not a path when supposed to be Hamiltonian? 
		if not has_cycle(G_pT, b[0]):
			if kills_branch2(b):
				stats[1] = stats[1] + skipped_count()
			else:
				for W in recurse():
					yield W
		# Restore variables
		G_pT[b[0]].remove(b[1])
		G_pT[b[1]].remove(b[0])
		rem_E.append(b)

	total_unfoldings = spanning_tree_count(G_T_rem_E)
	for W in recurse():
		yield W

# Input: a polycube
# Output: the area of the polycube
def surface_area(P):
	area = 0
	for cell in P:
		for neighvec in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
			neigh = (cell[0] + neighvec[0], cell[1] + neighvec[1], cell[2] + neighvec[2])
			if not neigh in P:
				area = area + 1
	return area	

# Input: a non-negative integer. 
# Output: a generator for the free (unique up to rotation & translation) polycubes
#         with the specified number of cells.
def enumerate_polycubes(n, cur=None):
	if cur == None:
		cur = n
	if cur < 1:
		return
	if cur == 1:
		yield set([(0, 0, 0)])
		return
	if cur == 2:
		yield set([(0, 0, 0), (0, 0, 1)])
		return
	if cur == 3:
		yield set([(0, 0, 0), (0, 0, 1), (0, 0, 2)])
		yield set([(0, 1, 0), (1, 0, 0), (0, 0, 0)])
		return

	def zero(P):
		zP = set([])
		delta = (min([x for x, y, z in P]), min([y for x, y, z in P]), min([z for x, y, z, in P]))
		for c in P:
			zP.add((c[0] - delta[0], c[1] - delta[1], c[2] - delta[2]))
		return zP

	def normalized_transformations(P):
		def rot(P, ci, cj, rep):
			rT = set([])
			for c in P:
				rc = c
				for i in xrange(rep):
					if set([ci, cj]) == set([0, 1]): 
						rc = (-rc[1], rc[0], rc[2])
					elif set([ci, cj]) == set([0, 2]): 
						rc = (-rc[2], rc[1], rc[0])
					else: 
						rc = (rc[0], -rc[2], rc[1])
				rT.add(rc)
			return rT
		
		trans = []
		for rep1 in xrange(4):
			for rep2 in xrange(4):
				trans.append(zero(rot(rot(P, 0, 1, rep1), 0, 2, rep2)))
			for rep3 in [1, 3]:	
				trans.append(zero(rot(rot(P, 0, 1, rep1), 1, 2, rep3)))
		return trans	

	already_enumerated = []
	for P in enumerate_polycubes(n, cur-1):
		for x in xrange(min([cx for cx, cy, cz in P])-1, max([cx for cx, cy, cz in P])+2):
			for y in xrange(min([cy for cx, cy, cz in P])-1, max([cy for cx, cy, cz in P])+2):
				for z in xrange(min([cz for cx, cy, cz in P])-1, max([cz for cx, cy, cz in P])+2):
					if (x, y, z) in P:
						continue
					cand = copy.copy(P)
					cand.add((x, y, z))
					if not is_polycube(cand):
						continue
					cand = zero(cand)
					found = False
					for aeP in already_enumerated:
						if cand == aeP:
							found = True
							break
					if not found:
						yield cand
						NT = normalized_transformations(cand)
						already_enumerated.extend(NT)


def enumerate_polycubes_with_area(area):
	for n in xrange(1, int(math.ceil((area-2)/4)) + 1):
		for P in enumerate_polycubes(n):
			if surface_area(P) == area:
				yield P

class TestStuff(unittest.TestCase):
	
	def setUp(self):
		pass

	def test__surface_area(self):
		self.assertEqual(6, surface_area(set([(0, 0, 0)])))
		self.assertEqual(10, surface_area(set([(0, 0, 0), (0, 0, 1)])))
		self.assertEqual(16, surface_area(set([(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1)])))

	def test__is_polycube(self):
		self.assertTrue(is_polycube(set([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])))
		self.assertTrue(is_polycube(set([(1, 1, 1), (1, 1, 2), (1, 2, 2)])))
		self.assertFalse(is_polycube(set([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)])))
		self.assertTrue(is_polycube(set([(0, 1, 0), (0, 2, 0), (1, 1, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)])))

	def test__unfoldings(self):
		# All unfoldings of cube tile
		count = 0
		for uf in unfoldings(set([(0, 0, 0)])):
			count = count + 1
			self.assertTrue(polyomino.is_polyomino(uf))
			self.assertTrue(isohedral.has_translation_tiling(uf) or isohedral.has_half_turn_tiling(uf))
		self.assertEqual(count, 384) # 384 = ((2*1)**3 * (2*2)**3 * (2*3)**1) / 8 

		# Lean on fact that every Hamiltonian unfolding has maximum boundary 
		for uf in unfoldings(set([(0, 0, 0)]), False, True):
			self.assertTrue(len(uf) == 2 + 2 * surface_area(set([(0, 0, 0)])))

	def test__unfolding_count(self):
		self.assertEqual(unfolding_count(set([(0, 0, 0)])), ((2*1)**3 * (2*2)**3 * (2*3)**1) / 8)

	def test__face_dual(self):
		V, E = face_dual(set([(0, 0, 0)]))
		self.assertEqual(len(V), 6) 
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

	def test__enumerate_polycubes(self):
		# https://oeis.org/A000162
		self.assertEqual(len([P for P in enumerate_polycubes(1)]), 1)
		self.assertEqual(len([P for P in enumerate_polycubes(2)]), 1)
		self.assertEqual(len([P for P in enumerate_polycubes(3)]), 2)
		self.assertEqual(len([P for P in enumerate_polycubes(4)]), 8)
		self.assertEqual(len([P for P in enumerate_polycubes(5)]), 29)

	def test__enumerate_polycubes_with_area(self):
		self.assertEqual(len([P for P in enumerate_polycubes_with_area(9)]), 0)
		self.assertEqual(len([P for P in enumerate_polycubes_with_area(10)]), 1)
		self.assertEqual(len([P for P in enumerate_polycubes_with_area(12)]), 0)
		

	def test__is_unfolding(self):
		cube = set([(0, 0, 0)])
		self.assertTrue(is_unfolding(cube, ['N', 'N', 'W', 'N', 'E', 'N', 'E', 'S', 'E', 'S', 'W', 'S', 'S', 'W']))
		self.assertTrue(is_unfolding(cube, ['E', 'E', 'N', 'E', 'S', 'E', 'S', 'W', 'S', 'W', 'N', 'W', 'W', 'N']))
		self.assertTrue(is_unfolding(cube, ['S', 'S', 'E', 'S', 'W', 'S', 'W', 'N', 'W', 'N', 'E', 'N', 'N', 'E']))
		self.assertTrue(is_unfolding(cube, ['N', 'N', 'E', 'S', 'E', 'E', 'E', 'S', 'S', 'W', 'N', 'W', 'W', 'W']))
		self.assertFalse(is_unfolding(cube, ['N', 'E', 'E', 'E', 'E', 'E', 'E', 'S', 'W', 'W', 'W', 'W', 'W', 'W']))
		
		F_pentacube_unfolding = ['W', 'N', 'E', 'N', 'N', 'N', 'E', 'S', 'S', 'E', 'S', 'E', 
			'S', 'S', 'E', 'S', 'S', 'S', 'S', 'W', 'W', 'W', 'N', 'E', 'E', 'N', 'N', 'W', 'N', 'N', 'S', 
			'S', 'E', 'S', 'W', 'W', 'N', 'W', 'W', 'N', 'E', 'E', 'N', 'N', 'S', 'W']

		F_pentacube = set([(0, 0, 0), (0, 1, 0), (-1, 1, 0), (0, 2, 0), (1, 2, 0)])
		for rot in [0, 90, 180, 270]:
			self.assertTrue(is_unfolding(F_pentacube, [polyomino.rot[rot][d] for d in F_pentacube_unfolding]))
			self.assertTrue(is_unfolding(F_pentacube, [polyomino.refl[0][polyomino.cw[d]] for d in F_pentacube_unfolding]))

		I_pentacube = set([(i, 0, 0) for i in xrange(5)])
		self.assertFalse(is_unfolding(I_pentacube, ['W', 'N', 'E', 'N', 'N', 'N', 'E', 'S', 'S', 'E', 'S', 'E', 
			'S', 'S', 'E', 'S', 'S', 'S', 'S', 'W', 'W', 'W', 'N', 'E', 'E', 'N', 'N', 'W', 'N', 'N', 'S', 
			'S', 'E', 'S', 'W', 'W', 'N', 'W', 'W', 'N', 'E', 'E', 'N', 'N', 'S', 'W']))

		dali_octacube = set([(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (-1, 0, 2), (1, 0, 2), (0, -1, 2), (0, 1, 2)])
		self.assertTrue(is_unfolding(dali_octacube, ['N', 'N', 'E', 'S', 'N', 'N', 'N', 'W', 'N', 'E', 'E', 
			'W', 'W', 'N', 'E', 'N', 'E', 'S', 'N', 'E', 'N', 'N', 'N', 'N', 'E', 'S', 'S', 'S', 'E', 
			'W', 'N', 'N', 'E', 'S', 'E', 'S', 'E', 'S', 'W', 'S', 'S', 'W', 'S', 'S', 'W', 'W', 'N', 
			'E', 'N', 'N', 'E', 'W', 'N', 'E', 'W', 'S', 'S', 'S', 'W', 'N', 'N', 'S', 'S', 'W', 'S', 
			'S', 'S', 'S', 'W', 'W']))

if __name__ == '__main__':
	unittest.main()

