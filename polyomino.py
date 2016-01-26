
import unittest

A = {'N', 'E', 'S', 'W'}

# Input: a list P of elements of A
# Output: whether P describes a closed path traversed clockwise 
# (ends where it begins)
def is_closed(P):
	result = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}
	current = (0, 0)
	for i in xrange(len(P)):
		current = (current[0] + result[P[i]][0], current[1] + result[P[i]][1])
	return current == (0, 0) 

# Input: a list P of elements of A
# Output: whether P describes a non-self-intersecting path
def is_simple(P):
	visited_locations = set([(0, 0)])
	result = {'N': (0, 1), 'E': (1, 0), 'S': (0, -1), 'W': (-1, 0)}
	current = (0, 0)
	for i in xrange(len(P)):
		current = (current[0] + result[P[i]][0], current[1] + result[P[i]][1])
		if current in visited_locations:
			# If not in special case of closed path 
			# finishing at start vertex
			if not (i == len(P)-1 and current == (0, 0)):
				return False
		visited_locations.add(current) 
	return True
	
# Input: a list P of elements of A
# Output: whether P is the boundary word of a polyomino
# (closed, simple, clockwise traversed)
def is_polyomino(P):
	if not is_closed(P) or not is_simple(P):
		return False
	result = {('N', 'E'): 1, ('N', 'N'): 0, ('N', 'W'): -1, 
		('E', 'S'): 1, ('E', 'E'): 0, ('E', 'N'): -1, 
		('S', 'W'): 1, ('S', 'S'): 0, ('S', 'E'): -1, 
		('W', 'N'): 1, ('W', 'W'): 0, ('W', 'S'): -1}
	winding_quart = result[(P[-1], P[0])] # Initialize with last turn
	for i in xrange(0, len(P)-1):
		winding_quart = winding_quart + result[(P[i], P[i+1])]
		current_dir = P[i]
	assert(winding_quart % 4 == 0) 
	return (winding_quart / 4 == 1)


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

	def test__is_polyomino(self):
		# Is ok
		P = ['N', 'E', 'S', 'W']
		for i in xrange(len(P)):
			self.assertTrue(is_polyomino(P[i:] + P[:i]))
		self.assertTrue(is_polyomino(['N', 'N', 'E', 'S', 'S', 'W']))
		self.assertTrue(is_polyomino(['N', 'N', 'W', 'N', 'E', 'E', 'S', 'S', 'S', 'S', 'W', 'N']))
		# Isn't clockwise
		P = ['N', 'W', 'S', 'E']
		for i in xrange(len(P)):
			self.assertFalse(is_polyomino(P[i:] + P[:i]))
		# Isn't simple		
		P = ['N', 'N', 'E', 'S', 'W', 'W', 'S', 'E']	
		for i in xrange(len(P)):	
			self.assertFalse(is_polyomino(P[i:] + P[:i])) 
		# Isn't closed
		self.assertFalse(is_polyomino(['N', 'E', 'S']))
		self.assertFalse(is_polyomino(['N', 'N', 'E', 'E', 'S', 'S', 'S', 'S', 'W', 'N'])) 
	
if __name__ == '__main__':
	unittest.main()

