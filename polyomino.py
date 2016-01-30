
import unittest
import copy

# The alphabet for the boundary words
A = {'N', 'E', 'S', 'W'}

# Some useful mappings
dir2vec = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}
comp = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}
vec2dir = {(0, 1): 'N', (1, 0): 'E', (0, -1): 'S', (-1, 0): 'W'}

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
	
# Input: a list of 2-tuples of numbers describing the vertices of a simple polygon 
# Output: the area of the polygon
def area(W):
	P = boundary_word2vertices(W)
	return 0.5 * abs(sum(P[i-1][0]*P[i][1] - P[i][0]*P[i-1][1] for i in xrange(len(P))))

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

class TestStuff(unittest.TestCase):

	def setUp(self):
		pass

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

	def test__area(self):
		self.assertEqual(area(['N', 'E', 'S', 'W']), 1)
		self.assertEqual(area(['N', 'N', 'E', 'S', 'S', 'W']), 2)
		self.assertEqual(area(['N', 'W', 'N', 'E', 'E', 'S', 'S', 'W']), 3)

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
		self.assertEqual(cancel(['N', 'E', 'N', 'S', 'S', 'W']), ['N', 'E', 'S', 'W'])
		self.assertEqual(cancel(['S', 'S', 'W', 'N', 'E', 'N']), ['S', 'W', 'N', 'E'])
		self.assertEqual(cancel(['N', 'E', 'N', 'S', 'S', 'E', 'W', 'W']), ['N', 'E', 'S', 'W'])
		W = ['N', 'S']
		cancel(W)
		self.assertEqual(W, ['N', 'S'])

if __name__ == '__main__':
	unittest.main()


