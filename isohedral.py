
import unittest
import polyomino


# Input: two words and a length.
# Output: the longest common prefix of the words, up to the specified length.
# Note: helper method
def longest_match(S1, S2, ub):
	i = 0
	while i < min(len(S1), len(S2), ub) and S1[i] == S2[i]:
		i = i + 1
	return i 	


# Input: a polyomino boundary word
# Output: a list of all admissible mirror factors
# Note: helper method
def admissible_mirror_factors(P):
	comp = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}

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
				factors.append(((i-l+n)%n, (i+r)%n))
	return factors


# Input: a polyomino boundary word
# Output: a list of all admissible mirror factors
# Note: helper method
def admissible_gapped_mirror_factor_pairs(P):
	comp = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E'}

	def inv_comp(S):
		return map(lambda l: comp[l], S[::-1])

	n = len(P)
	factor_pairs = []
	# Compute admissible mirror factors starting between letter pairs 
	for i in xrange(n):
		for j in xrange(i+1, n):
			l = longest_match(inv_comp(P + P[:i]), P[j:] + P, (i+n-j)/2)
			r = longest_match(P[i:] + P, inv_comp(P + P[:j]), (j-i)/2)
			if l == r and r > 0:
				start = (i-l+n) % n
				end = (i-1+r+n) % n
				factor_pairs.append((((i-l+n)%n, (i-1+r+n)%n), ((j-l+n)%n, (j-1+r+n)%n)))
	# Compute admissible mirror factors starting in middle of a letter
	for i in xrange(n):
		for j in xrange(i+1, n):
			if P[i] == comp[P[j]]:
				l = longest_match(inv_comp(P+P[:i]), P[(j+1)%n:] + P, (i+n-j-1)/2)
				r = longest_match(P[(i+1)%n:] + P, inv_comp(P + P[:j]), (j-i-1)/2) 
				if l == r:
					factor_pairs.append((((i-l+n)%n, (i+r)%n), ((j-l+n)%n, (j+r)%n)))
	return factor_pairs


# Input: a polyomino boundary word
# Output: a list of all admissible palindrome factors
# Note: helper method
def admissible_palindrome_factors(P):
	n = len(P)
	factors = []
	# Compute admissible palindrome factors starting between letter pairs 
	for i in xrange(n):
		l = longest_match((P + P[:i])[::-1], P[i:] + P, n/2)
		if l > 0:
			factors.append(((i-l+n)%n, (i-1+l)%n))
	# Compute admissible palindrome factors starting in middle of a letter
	for i in xrange(n):
		l = longest_match((P + P[:i])[::-1], P[(i+1)%n:] + P, n/2)
		factors.append(((i-l+n)%n, (i+l)%n))
	return factors

# Input: a polyomino boundary word
# Output: whether the polyomino has a quarter-turn tiling
def has_quarter_turn_tiling(P):
	pass		

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
				if AB_len == n/2:
					return True
				if ((B[1]+1)%n, (A[0]+n/2-1+n)%n) in factor_starts[(B[1]+1)%n]:
					return True
	return False

# Input: a polyomino boundary word
# Output: whether the polyomino has a translation tiling
def has_half_turn_tiling(P):
	n = len(P)
	mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P)
	palindrome_factor_start_lengths = {}
	palindrome_factor_end_lengths = {}
	for i in xrange(n):
		palindrome_factor_start_lengths[i] = []
		palindrome_factor_end_lengths[i] = []
	for f in admissible_palindrome_factors(P):
		f_length = f[1] + n * (f[1] < f[0]) - f[0] + 1
		palindrome_factor_start_lengths[f[0]].append(f_length)
		palindrome_factor_end_lengths[f[1]].append(f_length)	
	zero_length_mirror_factor_pairs = []
	for i in xrange(n):
		for j in xrange(i+1, n):
			zero_length_mirror_factor_pairs.append((((i+1)%n, i), ((j+1)%n, j)))
	for A, A_hat in mirror_factor_pairs + zero_length_mirror_factor_pairs:
		# Compute leftover double palindrome intervals
		dpi1 = ((A[1]+1)%n, (A_hat[0]-1+n)%n)
		dpi2 = ((A_hat[1]+1)%n, (A[0]-1+n)%n)
		# Check if both intervals are double palindromes
		dpis = [dpi1, dpi2]
		dpi_pass = {dpi1: False, dpi2: False}
		for dpi in dpis:
			for B_len in palindrome_factor_start_lengths[dpi[0]]:
				dpi_len = dpi[1] + n * (dpi[1] < dpi[0]) - dpi[0] + 1
				if dpi_len == B_len:
					dpi_pass[dpi] = True
					break
				for C_len in palindrome_factor_end_lengths[dpi[1]]:
					if dpi_len == B_len + C_len:
						dpi_pass[dpi] = True
						break
		if dpi_pass[dpi1] and dpi_pass[dpi2]:
			return True	
	return False

class TestStuff(unittest.TestCase):

	def setUp(self):
		pass
	
	def test__admissible_palindrome_factors(self):
		self.assertEqual(set(admissible_palindrome_factors(['N', 'E', 'S', 'W'])), set([(0, 0), (1, 1), (2, 2), (3, 3)])) 

	def test__has_translation_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 1xi, ix1 rectangles
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N'] + ['E']*i + ['S'] + ['W']*i))
		# Small tetris pieces (pseudosquare tiling)
		self.assertTrue(has_translation_tiling(['N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'N', 'W']))
		# Small tetris pieces (pseudohex tiling)
		self.assertTrue(has_translation_tiling(['N', 'N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'W']))
		self.assertTrue(has_translation_tiling(['N', 'E']*3 + ['S', 'E', 'E'] + ['S', 'W']*2 + 
			['N', 'W', 'S', 'W', 'W']))
		self.assertFalse(has_translation_tiling(['N', 'E']*3 + ['S', 'E', 'E'] + ['S', 'S', 'W', 'W'] + 
			['N', 'W', 'S', 'W', 'W']))
		# Small tetris pieces (no tiling)
		self.assertFalse(has_translation_tiling(['N', 'E', 'N', 'N', 'E', 'S', 'S', 'S', 
			'E', 'E', 'S', 'W', 'W', 'W', 'N', 'W']))
		self.assertFalse(has_translation_tiling(['W', 'W', 'W', 'N', 'N', 'E', 'S', 'E', 'N', 'E', 'S', 'S']))	
	
	def test__has_half_turn_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_half_turn_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 1xi, ix1 rectangles
		for i in [1,2,3,4,5]:
			self.assertTrue(has_half_turn_tiling(['N'] + ['E']*i + ['S'] + ['W']*i))
		# Small tetris pieces (one palindrome on each side)
		self.assertTrue(has_half_turn_tiling(['N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'N', 'W']))
		self.assertTrue(has_half_turn_tiling(['N', 'N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'W']))
		# Small tetris pieces (one palindrome on one side, two on the other)
		self.assertTrue(has_half_turn_tiling(['N', 'W', 'N', 'E', 'E', 'N', 'W', 'N'] +
			['E', 'N', 'E', 'S', 'E'] +
			['S', 'E', 'S'] +
			['S', 'W', 'S'] +
			['W', 'N', 'W', 'S', 'W']))
		# Small tetris pieces (two palindromes on both sides)
		self.assertTrue(has_half_turn_tiling(['N', 'E']*3 + ['S', 'E', 'E'] + ['S', 'W']*2 + ['N', 'W', 'S', 'W', 'W']))
		self.assertTrue(has_half_turn_tiling(['W', 'W', 'W', 'N', 'N', 'E', 'S', 'E', 'N', 'E', 'S', 'S']))	
		# Small tetris pieces (no translation section)
		B = ['E', 'N', 'E', 'S', 'E', 'E', 'N', 'E', 'S']
		B = B + B[::-1]
		C = ['S', 'S', 'E', 'S', 'S', 'W']
		C = C + C[::-1]
		D = ['W', 'W', 'W', 'W', 'W']
		D = D + D[::-1]
		E = ['N', 'W', 'W', 'N', 'E', 'E']
		E = E + E[::-1]
		self.assertTrue(has_half_turn_tiling(B + C + D + E))
		# Small tetris pieces (no tilings)
		self.assertFalse(has_half_turn_tiling(['N', 'E', 'N', 'N', 'E', 'S', 'S', 'S', 
			'E', 'E', 'S', 'W', 'W', 'W', 'N', 'W']))
		self.assertFalse(has_half_turn_tiling(['W']*5 + ['N', 'N', 'E', 'S'] + ['E']*3 + ['N', 'E', 'S', 'S']))

if __name__ == '__main__':
	unittest.main()



