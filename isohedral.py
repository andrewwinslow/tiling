
import unittest
import polyomino

# Input: a polyomino boundary word
# Output: a list of all admissible mirror factors
# Note: helper method
def admissible_mirror_factors(P):
	comp = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}

	def longest_match(S1, S2, ub):
		i = 0
		while i < min(len(S1), len(S2), ub) and S1[i] == S2[i]:
			i = i + 1
		return i 	
	
	def inv_comp(S):
		return map(lambda l: comp[l], S[::-1])

	n = len(P)
	factors = []
	# Compute admissible mirror factors starting between letter pairs 
	for i in xrange(n):
		l = longest_match(inv_comp(P + P[:i]), P[(i+n/2)%n:] + P, n/4)
		r = longest_match(P[i:] + P, inv_comp(P + P[:(i+n/2)%n]), n/4)
		if l == r and r > 0:
			start = (i-l+n) % n
			end = (i-1+r+n) % n
			factors.append((start, end))
	# Compute admissible mirror factors starting in middle of a letter
	for i in xrange(n):
		if P[i] == comp[P[(i+n/2) % n]]:
			l = longest_match(inv_comp(P+P[:i]), P[(i+n/2+1)%n:] + P, (n-2)/4)
			r = longest_match(P[(i+1)%n:] + P, inv_comp(P + P[:(i+n/2)%n]), (n-2)/4) 
			if l == r:
				start = (i-l+n) % n
				end = (i+r) % n
				factors.append((start, end))
	return factors


# Input: a polyomino boundary word
# Output: whether the polyomino has a translation tiling
def has_translation_tiling(P):
	n = len(P)
	factors = admissible_mirror_factors(P)
	factor_starts = {}
	factor_ends = {}
	for i in xrange(n):
		factor_starts[i] = set([])
		factor_ends[i] = set([])
	for f in factors:
		factor_starts[f[0]].add(f)
		factor_ends[f[1]].add(f)
	for i in xrange(n):
		for A in factor_starts[i]:
			for B in factor_starts[(A[1]+1)%n]:
				AB_len = (A[1]-A[0]+1+n) % n + (B[1]-B[0]+1+n) % n 
				if AB_len > n/2:
					continue
				if AB_len == n/2 or ((B[1]+1)%n, (A[0]-1+n)%n) in factor_starts[(B[1]+1)%n]:
					return True


class TestStuff(unittest.TestCase):

	def setUp(self):
		pass

	def test__has_translation_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 1xi, ix1 rectangles
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N'] + ['E']*i + ['S'] + ['W']*i))
		# Small tetris piece 1 (pseudosquare tiling)
		self.assertTrue(has_translation_tiling(['N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'N', 'W']))
		self.assertFalse(has_translation_tiling(['N', 'N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'W']))
		# Small tetris piece 2 (pseudohex tiling)
		self.assertTrue(['N', 'E']*3 + ['S', 'E', 'E'] + ['S', 'W']*2 + ['N', 'W', 'S', 'W', 'W'])
		self.assertTrue(['N', 'E']*3 + ['S', 'E', 'E'] + ['S', 'S', 'W', 'W'] + ['N', 'W', 'S', 'W', 'W'])
		
	def test__has_half_turn_tiling(self):
		pass

if __name__ == '__main__':
	unittest.main()

