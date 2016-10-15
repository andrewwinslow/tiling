
import unittest
import copy
import random 

# The alphabet for the boundary words
A = set(['N', 'E', 'S', 'W'])

# Some useful mappings
dir2vec = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}
comp = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}
vec2dir = {(0, 1): 'N', (1, 0): 'E', (0, -1): 'S', (-1, 0): 'W'}
cw = {'N': 'E', 'E': 'S', 'S': 'W', 'W': 'N'}
ccw = {'N': 'W', 'E': 'N', 'S': 'E', 'W': 'S'}
rot = {0: {'N': 'N', 'E': 'E', 'S': 'S', 'W': 'W'},
	90: {'N': 'W', 'W': 'S', 'S': 'E', 'E': 'N'},
	180: {'N': 'S', 'W': 'E', 'S': 'N', 'E': 'W'},
	270: {'N': 'E', 'W': 'N', 'S': 'W', 'E': 'S'}} 
refl = {-45: {'N': 'W', 'E': 'S', 'S': 'E', 'W': 'N'}, 
	0: {'N': 'S', 'E': 'E', 'S': 'N', 'W': 'W'},
	45: {'N': 'E', 'E': 'N', 'S': 'W', 'W': 'S'},
	90: {'N': 'N', 'E': 'W', 'S': 'S', 'W': 'E'}} 



# Input: a list W of element of A
# Output: W with adjacent opposite direction elements removed
def cancel(W):
	W = copy.copy(W)
	cancel = True
	while cancel and len(W) >= 2:
		cancel = False
		for i in xrange(-1, len(W)-1):
			if W[i] == comp[W[i+1]]:
				if i != -1:
					del W[i:i+2]
				else:
					del W[-1]
					del W[0]
				cancel = True
				break
	return W

# Input: a list W of elements of A
# Output: whether W describes a closed path traversed clockwise 
# (ends where it begins)
def is_closed(W):
	current = (0, 0)
	for i in xrange(len(W)):
		current = (current[0] + dir2vec[W[i]][0], current[1] + dir2vec[W[i]][1])
	return current == (0, 0) 

# Input: a polyomino boundary word
# Output: the polyomino as a set of cells specified as integer 2-tuples
def word2polyomino(W):
	column2edges = {}
	loc = (0, 0)
	for i in xrange(len(W)):
		if W[i] == 'E':
			if not loc[0] in column2edges:
				column2edges[loc[0]] = [] 
			column2edges[loc[0]].append(loc[1])
		if W[i] == 'W':
			if not loc[0]-1 in column2edges:
				column2edges[loc[0]-1] = [] 
			column2edges[loc[0]-1].append(loc[1])
		loc = (loc[0] + dir2vec[W[i]][0], loc[1] + dir2vec[W[i]][1])
	cells = []
	for x in column2edges:
		for y in xrange(min(column2edges[x]), max(column2edges[x])+1):
			if len(filter(lambda e: e > y, column2edges[x])) % 2 == 1:
				cells.append((x, y))
	return set(cells)	 

# Input: a polyomino boundary word W
# Output: whether W describes an orthogonally convex polyomino
def is_orthoconvex(W):
	verts = boundary_word2vertices(W)
	x_coords = [v[0] for v in verts]
	y_coords = [v[1] for v in verts]

	i = -1
	while i < len(verts) and x_coords[i-1] <= x_coords[i]: 
		i = i + 1
	while i < len(verts) and x_coords[i-1] >= x_coords[i]: 
		i = i + 1
	while i < len(verts) and x_coords[i-1] <= x_coords[i]: 
		i = i + 1
	if i != len(verts):
		return False

	i = -1
	while i < len(verts) and y_coords[i-1] <= y_coords[i]: 
		i = i + 1
	while i < len(verts) and y_coords[i-1] >= y_coords[i]: 
		i = i + 1
	while i < len(verts) and y_coords[i-1] <= y_coords[i]: 
		i = i + 1
	if i != len(verts):
		return False
	
	return True	

# Input: a list W of elements of A
# Output: whether W describes a non-self-intersecting path
def is_simple(W):
	# Pass #1: any degenerate opposite directions following each other?
	for i in xrange(len(W)-1):
		if W[i+1] == comp[W[i]]:
			return False
	# Pass #2: any points revisited (except (0, 0) as last vertex)
	visited_locations = set([(0, 0)])
	current = (0, 0)
	for i in xrange(len(W)):
		current = (current[0] + dir2vec[W[i]][0], current[1] + dir2vec[W[i]][1])
		if current in visited_locations:
			# If not in special case of closed path 
			# finishing at start vertex
			if not (i == len(W)-1 and current == (0, 0)):
				return False
		visited_locations.add(current) 
	return True

# Input: a closed simple boundary word
# Output: whether the boundary word is traversed clockwise
# Note: helper method
def is_clockwise(W):
	if len(W) == 0:
		return True
	result = {('N', 'E'): 1, ('N', 'N'): 0, ('N', 'W'): -1, 
		('E', 'S'): 1, ('E', 'E'): 0, ('E', 'N'): -1, 
		('S', 'W'): 1, ('S', 'S'): 0, ('S', 'E'): -1, 
		('W', 'N'): 1, ('W', 'W'): 0, ('W', 'S'): -1}
	winding_quart = result[(W[-1], W[0])] # Initialize with last turn
	for i in xrange(0, len(W)-1):
		winding_quart = winding_quart + result[(W[i], W[i+1])]
		current_dir = W[i]
	return (winding_quart / 4 == 1)

# Input: a list W of elements of A
# Output: whether W is the boundary word of a polyomino
# (closed, simple, clockwise traversed)
def is_polyomino(W):
	return is_closed(W) and is_simple(W) and is_clockwise(W)

# Input: a polyomino boundary word
# Output: the vertices of a polyomino with this boundary word
def boundary_word2vertices(W):
	current = (0, 0)
	verts = [(0, 0)]
	for i in xrange(len(W)-1):
		current = (current[0] + dir2vec[W[i]][0], current[1] + dir2vec[W[i]][1])
		verts.append(current)
	return verts 

# Input: a lattice path
# Output: the vertices of the lattice path
def lattice_path2vertices(W):
	current = (0, 0)
	verts = [(0, 0)]
	for i in xrange(len(W)):
		current = (current[0] + dir2vec[W[i]][0], current[1] + dir2vec[W[i]][1])
		verts.append(current)
	return verts 
	
# Input: a list of 2-tuples of numbers describing the vertices of a simple polygon 
# Output: the area of the polygon
def area(W):
	P = boundary_word2vertices(W)
	return int(0.5 * abs(sum(P[i-1][0]*P[i][1] - P[i][0]*P[i-1][1] for i in xrange(len(P)))))

# Input: a polyomino boundary word
# Output: whether the polyomino is weakly simple (shares vertices but no overlapping cells)
def is_weakly_simple(W):
	cells = word2polyomino(W)
	return len(cells) == area(W)	

# Input: a positive integer n
# Output: generator of the boundary words of length n
def enumerate_boundary_words(n):

	def neighbors(p):
		return [(p[0], p[1]+1), (p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]-1)]

	def dir(p1, p2):
		return (p2[0]-p1[0], p2[1]-p1[1])

	path = [(0, 0)]
	def recurse():
		if len(path) > n:
			return
		if len(path) == n:
			if not ((0, 0) in neighbors(path[-1])):
				return 
			# Convert to a boundary word and yield it
			word = (map(lambda i: vec2dir[dir(path[i], path[i+1])], xrange(len(path)-1)) + 
				[vec2dir[dir(path[-1], (0, 0))]])
			if is_polyomino(word) and min(path) == (0, 0):
				yield word
		head = path[-1]
		for new in [(head[0], head[1]+1), (head[0]+1, head[1]), (head[0]-1, head[1]), (head[0], head[1]-1)]:
			if not (new in path): # Check for self-intersection
				path.append(new)
				for w in recurse():
					yield w
				path.pop() 
	for w in recurse():
		yield w

# Input: a length l
# Output: all strongly simple paths of length l
def enumerate_strongly_simple_paths(l):
	if l == 0:
		yield []
		return
	W = []
	def recurse(W):
		if len(W) == l:
			yield W
			return
		for next_d in A:
			if next_d == comp[W[-1]]:
				continue
			W.append(next_d)
			if len(set(lattice_path2vertices(W))) == len(W)+1:
				for RW in recurse(W):
					yield RW	
			W.pop()
	for s in A:
		W.append(s)
		for RW in recurse(W):
			yield RW
		W.pop()

# Input: a length l
# Output: all strongly simple palindromic paths of length l
def enumerate_strongly_simple_palindromes(l):	
	if l == 0:
		yield []
		return
	if l == 1:
		for d in A:
			yield [d]
		return
	W = []
	def recurse(W):
		if not len(W) < l/2:
			if l % 2 == 0:
				yield W[::-1] + W
			else: # Odd length case
				for next_d in A:
					if next_d == comp[W[0]]:
						continue
					cand = W[::-1] + [next_d] + W
					if len(set(lattice_path2vertices(cand))) == len(cand)+1:
						yield cand 
			return
		for next_d in A:
			if next_d == comp[W[-1]]:
				continue
			W.append(next_d)
			if len(set(lattice_path2vertices(W[::-1] + W))) == 2*len(W)+1: 
				for RW in recurse(W):
					yield RW
			W.pop() 

	for s in A:
		W.append(s)
		for RW in recurse(W):
			yield RW
		W.pop()

# Input: a length l
# Output: all strongly simple double palindromic paths of length l
def enumerate_strongly_simple_double_palindromes(l):
	for l1 in xrange(l+1):
		for l2 in xrange(l-1-l1):
			for p1 in enumerate_strongly_simple_palindromes(l1):
				for p2 in enumerate_strongly_simple_palindromes(l2):
					cand = p1 + p2
					if len(set(lattice_path2vertices(cand))) == len(cand)+1:
						yield cand 
	for p in enumerate_strongly_simple_palindromes(l):
		yield p

# Input: a length l. Optionally, an area a.
# Output: all strongly simple boundary words of length l (and area a) obeying conways criterion.
def enumerate_conways_boundary_words(l, a=None):
	vec2dp = {}
	vecs = set([])
	for dp in enumerate_strongly_simple_double_palindromes(l):
		verts = lattice_path2vertices(dp)
		vec = (verts[-1][0] - verts[0][0], verts[-1][1] -  verts[0][1])
		if not (vec, len(dp)) in vec2dp:
			vec2dp[(vec, len(dp))] = []
		vec2dp[(vec, len(dp))].append(dp)
		vecs.add(vec)

	def enumerate_partitions(rem, p, l23min):
		if len(p) == 0:
			for i in xrange(2*rem+1):
				p.append(i)
				for r in enumerate_partitions(rem-2*i, p, l23min):
					yield r
				p.pop()
		elif len(p) < 3: # l2, l3	
			for i in xrange(l23min, rem+1):
				p.append(i)
				for r in enumerate_partitions(rem-i, p, l23min):
					yield r
				p.pop()
		else:
			if p[0]*2 + p[1] + p[2] == l:
				yield p	
	
	def mag(vec):
		return abs(vec[0]) + abs(vec[1])

	for vec in vecs:
		nvec = (-vec[0], -vec[1])
		for part in enumerate_partitions(l, [], mag(vec)):
			l1 = part[0]
			l2 = part[1]
			l3 = part[2]
			if (not (vec, l2) in vec2dp) or (not (nvec, l3) in vec2dp):
				continue 
			for A in enumerate_strongly_simple_paths(l1):
				A_hat = [comp[s] for s in A[::-1]]
				for dp1 in vec2dp[(vec, l2)]:
					if not is_simple(A + dp1 + A_hat):
						continue
					for dp2 in vec2dp[(nvec, l3)]:
						cand = A + dp1 + A_hat + dp2
						if area != None and area(cand) != a:
							continue
						if is_polyomino(cand):
							yield cand

# Input: a weakly simple polyomino boundary word
# Output: the cell dual	graph as a vertex set and edge set	
def cell_dual(W):
	column2edges = {}
	row2edges = {}
	loc = (0, 0)
	for i in xrange(len(W)):
		if W[i] == 'E':
			if not loc[0] in column2edges:
				column2edges[loc[0]] = [] 
			column2edges[loc[0]].append(loc[1])
		elif W[i] == 'W':
			if not loc[0]-1 in column2edges:
				column2edges[loc[0]-1] = [] 
			column2edges[loc[0]-1].append(loc[1])
		elif W[i] == 'N':
			if not loc[1] in row2edges:
				row2edges[loc[1]] = []
			row2edges[loc[1]].append(loc[0])
		else: # W[i] == 'S':
			if not loc[1]-1 in row2edges:
				row2edges[loc[1]-1] = []
			row2edges[loc[1]-1].append(loc[0])
		loc = (loc[0] + dir2vec[W[i]][0], loc[1] + dir2vec[W[i]][1])
	V = set([])
	for x in column2edges:
		for y in xrange(min(column2edges[x]), max(column2edges[x])+1):
			if len(filter(lambda e: e > y, column2edges[x])) % 2 == 1:
				V.add((x, y))
	E = set([])
	for c in V:
		neigh = (c[0] + 0, c[1] + 1)
		if neigh in V and not (neigh[1] in column2edges[neigh[0]]):
			E.add((c, neigh))
		neigh = (c[0] + 1, c[1])
		if neigh in V and not (neigh[0] in row2edges[neigh[1]]):
			E.add((c, neigh))
	return V, E

class TestStuff(unittest.TestCase):

	def setUp(self):
		pass

	def test__cell_dual(self):
		self.assertEqual(cell_dual(['N', 'N', 'N', 'E', 'E', 'E', 'S', 'S', 'S', 'W', 'W', 'W']),
			(set([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]), 
				set([((0, 0), (1, 0)), ((0, 0), (0, 1)), ((1, 0), (2, 0)), ((1, 0), (1, 1)), ((2, 0), (2, 1)),
					((0, 1), (1, 1)), ((0, 1), (0, 2)), ((1, 1), (2, 1)), ((1, 1), (1, 2)),
					((2, 1), (2, 2)), ((0, 2), (1, 2)), ((1, 2), (2, 2))])))
		self.assertEqual(cell_dual(['N', 'E', 'W', 'N', 'E', 'E', 'S', 'S', 'W', 'W']), 
			(set([(0, 0), (0, 1), (1, 0), (1, 1)]), set([((0, 0), (1, 0)), ((1, 0), (1, 1)), ((0, 1), (1, 1))])))

	def test__is_closed(self):
		# Clockwise
		self.assertTrue(is_closed(['N', 'E', 'S', 'W']))
		self.assertTrue(is_closed(['S', 'W', 'N', 'E']))
		self.assertTrue(is_closed(['N', 'N', 'E', 'S', 'W', 'W', 'S', 'E']))
		# Counterclockwise
		self.assertTrue(is_closed(['N', 'W', 'S', 'E']))
		self.assertTrue(is_closed(['S', 'E', 'N', 'W']))
		self.assertTrue(is_closed(['N', 'N', 'W', 'S', 'E', 'E', 'S', 'W']))
		# Just don't make it
		self.assertFalse(is_closed(['N', 'E', 'S', 'S']))	
		# Degenerate
		self.assertTrue(is_closed(['N', 'N', 'S', 'S']))

	def test__is_simple(self):
		# Simple but possible closed
		self.assertTrue(is_simple(['N', 'E', 'S']))
		self.assertTrue(is_simple(['N', 'E', 'S', 'W']))
		self.assertTrue(is_simple(['N', 'W', 'S', 'E']))
		self.assertTrue(is_simple(['S', 'E', 'N', 'W']))
		# Degenerate
		self.assertFalse(is_simple(['N', 'S']))
		self.assertFalse(is_simple(['N', 'E', 'S', 'N']))
		self.assertFalse(is_simple(['N', 'N', 'N', 'S']))
		# Typical self-intersection
		self.assertFalse(is_simple(['N', 'N', 'E', 'S', 'W', 'W']))
		self.assertFalse(is_simple(['N', 'N', 'W', 'S', 'E', 'E']))
		# Self-intersection because it's closed plus some edges
		self.assertFalse(is_simple(['N', 'E', 'S', 'W', 'W']))
		self.assertFalse(is_simple(['N', 'E', 'S', 'W', 'N']))
		self.assertFalse(is_simple(['N', 'W', 'S', 'E', 'E']))
		self.assertFalse(is_simple(['N', 'W', 'S', 'E', 'N']))
		# Self-intersection because it's closed but retraces the path
		self.assertFalse(is_simple(['N', 'E', 'S', 'W']*2))
		pass

	def test_is_clockwise(self):
		W = ['N', 'E', 'S', 'W']
		for i in xrange(len(W)):
			self.assertTrue(is_polyomino(W[i:] + W[:i]))
		W = ['N', 'W', 'S', 'E']
		for i in xrange(len(W)):
			self.assertFalse(is_polyomino(W[i:] + W[:i]))

	def test__is_polyomino(self):
		# Is ok
		W = ['N', 'E', 'S', 'W']
		for i in xrange(len(W)):
			self.assertTrue(is_polyomino(W[i:] + W[:i]))
		self.assertTrue(is_polyomino(['N', 'N', 'E', 'S', 'S', 'W']))
		self.assertTrue(is_polyomino(['N', 'N', 'W', 'N', 'E', 'E', 'S', 'S', 'S', 'S', 'W', 'N']))
		self.assertTrue(is_polyomino([]))
		# Isn't clockwise
		W = ['N', 'W', 'S', 'E']
		for i in xrange(len(W)):
			self.assertFalse(is_polyomino(W[i:] + W[:i]))
		# Isn't simple		
		W = ['N', 'N', 'E', 'S', 'W', 'W', 'S', 'E']	
		for i in xrange(len(W)):	
			self.assertFalse(is_polyomino(W[i:] + W[:i])) 
		# Isn't closed
		self.assertFalse(is_polyomino(['N', 'E', 'S']))
		self.assertFalse(is_polyomino(['N', 'N', 'E', 'E', 'S', 'S', 'S', 'S', 'W', 'N'])) 

	def test__boundary_word2vertices(self):
		self.assertEqual(boundary_word2vertices(['N', 'E', 'S', 'W']), [(0, 0), (0, 1), (1, 1), (1, 0)])
		self.assertEqual(boundary_word2vertices(['N', 'N', 'E', 'S', 'S', 'W']), 
			[(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0)])
		self.assertEqual(boundary_word2vertices(['N', 'E', 'E', 'S', 'W', 'W']), 
			[(0, 0), (0, 1), (1, 1), (2, 1), (2, 0), (1, 0)])
		self.assertEqual(boundary_word2vertices(
			['N', 'N', 'N', 'E', 'E', 'E', 'S', 'W', 'W', 'S', 'E', 'N', 'E', 'S', 'S', 'W', 'W', 'W']),
			[(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (3, 3), (3, 2), (2, 2), (1, 2), 
				(1, 1), (2, 1), (2, 2), (3, 2), (3, 1), (3, 0), (2, 0), (1, 0)])

	def test__area(self):
		self.assertEqual(area(['N', 'E', 'S', 'W']), 1)
		self.assertEqual(area(['N', 'N', 'E', 'S', 'S', 'W']), 2)
		self.assertEqual(area(['N', 'W', 'N', 'E', 'E', 'S', 'S', 'W']), 3)
		self.assertEqual(area(['N', 'N', 'N', 'E', 'E', 'E', 'S', 'W', 'W', 'S', 'E', 'N', 'E', 'S', 'S', 'W', 'W', 'W']), 8)

	def test__enumerate_boundary_words(self):
		for i in xrange(4):
			self.assertEqual(len([W for W in enumerate_boundary_words(i)]), 0)

		result = [W for W in enumerate_boundary_words(4)]
		self.assertEqual(len(result), 1)
		W = ['N', 'E', 'S', 'W']
		self.assertTrue(W in result)

		result = [W for W in enumerate_boundary_words(6)]
		self.assertEqual(len(result), 2)
		W = ['N', 'N', 'E', 'S', 'S', 'W']
		self.assertTrue(W in result)
		W = ['N', 'E', 'E', 'S', 'W', 'W']
		self.assertTrue(W in result)

	def test__cancel(self):
		self.assertEqual(cancel(['S', 'N']), [])
		self.assertEqual(cancel(['S', 'S', 'N', 'N']), [])
		self.assertEqual(cancel(['E', 'W', 'N', 'N', 'S', 'S', 'W', 'E', 'E', 'E']), ['E', 'E'])
		self.assertEqual(cancel(['N', 'E', 'N', 'S', 'S', 'W']), ['N', 'E', 'S', 'W'])
		self.assertEqual(cancel(['S', 'S', 'W', 'N', 'E', 'N']), ['S', 'W', 'N', 'E'])
		self.assertEqual(cancel(['N', 'E', 'N', 'S', 'S', 'E', 'W', 'W']), ['N', 'E', 'S', 'W'])
		self.assertEqual(cancel(['N', 'N', 'E', 'E', 'S', 'S']), ['E', 'E'])
		self.assertEqual(cancel(['N', 'N', 'W', 'S', 'E', 'E', 'S', 'S']), ['S', 'E'])
		for i in xrange(100):
			instance = [random.choice(list(A)) for j in xrange(500)]
			cancelled = cancel(instance)
			for j in xrange(-1, len(cancelled)-1):
				self.assertNotEqual(cancelled[j], comp[cancelled[j+1]])
		W = ['N', 'S']
		cancel(W)
		self.assertEqual(W, ['N', 'S'])
		self.assertEqual(cancel(['N', 'S', 'E', 'E', 'N', 'S', 'E', 'N', 'S']), ['E', 'E', 'E'])

	def test__word2polyomino(self):
		self.assertEqual(word2polyomino(['N', 'E', 'S', 'W']), set([(0, 0)]))
		self.assertEqual(word2polyomino(['N', 'N', 'E', 'S', 'S', 'W']), set([(0, 0), (0, 1)]))
		self.assertEqual(word2polyomino(['N', 'N', 'E', 'S', 'E', 'S', 'W', 'W']), set([(0, 0), (0, 1), (1, 0)]))
		self.assertEqual(word2polyomino(['N', 'N', 'E', 'S', 'E', 'N', 'E', 'S', 'S', 'W', 'W', 'W']), 
			set([(0, 0), (0, 1), (1, 0), (2, 0), (2, 1)]))
		self.assertEqual(word2polyomino(['N', 'N', 'E', 'E', 'S', 'W', 'S', 'E', 'S', 'W', 'W', 'N']),
			set([(0, 0), (0, 1), (1, 1), (0, -1), (1, -1)]))

	def test__is_weakly_simple(self):
		self.assertTrue(is_weakly_simple(['N', 'E', 'S', 'W']))
		self.assertTrue(is_weakly_simple(['N', 'N', 'E', 'S', 'S', 'W']))
		self.assertTrue(is_weakly_simple(['N', 'N', 'E', 'E', 'S', 'W', 'E', 'S', 'W', 'W']))
		self.assertTrue(is_weakly_simple(['N', 'N', 'N', 'E', 'E', 'E', 'S', 'W', 'W', 'S', 'E', 'N', 'E', 'S', 'S', 'W', 'W', 'W']))
		self.assertFalse(is_weakly_simple(['N', 'N', 'E', 'E', 'E', 'S', 'W', 'W', 'S', 'E', 'N', 'N', 'E', 'S', 'S', 'S', 'W', 'W', 'W']))

if __name__ == '__main__':
	unittest.main()


