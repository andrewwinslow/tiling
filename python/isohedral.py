
import unittest
import polyomino

# Input: two words and a length.
# Output: the longest common prefix of the words, up to the specified length.
def longest_match(S1, S2, ub):
	i = 0
	while i < min(len(S1), len(S2), ub) and S1[i] == S2[i]:
		i = i + 1
	return i 	

# Input: a polyomino boundary word
# Output: a list of all admissible mirror factors
def admissible_mirror_factors(P):

	def inv_comp(S):
		return map(lambda l: polyomino.comp[l], S[::-1])

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
		if P[i] == polyomino.comp[P[(i+n/2) % n]]:
			l = longest_match(inv_comp(P+P[:i]), P[(i+n/2+1)%n:] + P, (n-2)/4)
			r = longest_match(P[(i+1)%n:] + P, inv_comp(P + P[:(i+n/2)%n]), (n-2)/4) 
			if l == r:
				factors.append(((i-l+n)%n, (i+r)%n))
	return factors


# Input: a polyomino boundary word
# Output: a list of all admissible mirror factors
# Note: only gives one tuple for each pair 
#       (the one whose center has lower index) 
def admissible_gapped_mirror_factor_pairs(P):

	def inv_comp(S):
		return map(lambda l: polyomino.comp[l], S[::-1])

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
			if P[i] == polyomino.comp[P[j]]:
				l = longest_match(inv_comp(P+P[:i]), P[(j+1)%n:] + P, (i+n-j-1)/2)
				r = longest_match(P[(i+1)%n:] + P, inv_comp(P + P[:j]), (j-i-1)/2) 
				if l == r:
					factor_pairs.append((((i-l+n)%n, (i+r)%n), ((j-l+n)%n, (j+r)%n)))
	return factor_pairs

# Input: a polyomino boundary word and angle theta in set([-45, 0, 45, 90])
# Output: the admissible gapped reflect squares of angle theta
# Note: only gives one tuple for each pair 
#       (the one whose start has lower index)
def admissible_gapped_reflect_square_factor_pairs(P, theta):
	n = len(P)
	factor_pairs = []
	for i in xrange(n):
		for j in range(i+1, n):
			d = min(j - i + n * (j < i), i - j + n * (i < j))
			l = longest_match(P[i:] + P, [polyomino.refl[theta][s] for s in P[j:] + P], d+1)
			if 1 <= l <= d:
				factor_pairs.append(((i, (i+l-1+n)%n), (j, (j+l-1+n)%n)))
	return factor_pairs 

# Input: a polyomino boundary word 
# Output: the admissible reflect squares of the word (for all theta)
def admissible_reflect_square_factors(P):
	n = len(P)

	def is_reflect_square_factor(i, j, theta):
		l = j - i + 1 + n * (j < i)
		if l % 2 != 0:
			return False
		l = l / 2 # Change to repeat len
		if l == longest_match(P[i:] + P, [polyomino.refl[theta][s] for s in P[i+l:] + P], l+1):
			return True

	factors = []
	for i in xrange(n):
		for j in range(i) + range(i+1, n):
			for theta in [-45, 0, 45, 90]:
				if is_reflect_square_factor(i, j, theta):
					factors.append((i, j)) 
	return list(set(factors))	

# Input: a polyomino boundary word and angle theta
# Output: a list of all admissible theta-drome factors
def admissible_rotadrome_factors(P, theta):
	rot_f = {180: lambda x: x, 90: lambda x: polyomino.ccw[x]}
	n = len(P)
	factors = []
	# Compute admissible palindrome factors starting between letter pairs 
	for i in xrange(n):
		l = longest_match((P + P[:i])[::-1], [rot_f[theta](s) for s in P[i:] + P], n/2)
		if l > 0:
			factors.append(((i-l+n)%n, (i-1+l)%n))
	# Compute admissible palindrome factors starting in middle of a letter
	if theta == 180:
		for i in xrange(n):
			l = longest_match((P + P[:i])[::-1], P[(i+1)%n:] + P, n/2)
			factors.append(((i-l+n)%n, (i+l)%n))
	return factors

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_type_1_reflection_tiling(P):
	n = len(P)
	mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P)
	reflect_square_factors = admissible_reflect_square_factors(P)
	reflect_square_factor_starts = [[] for i in xrange(n)]
	for f in reflect_square_factors:
		reflect_square_factor_starts[f[0]].append(f)
	for A, A_hat in mirror_factor_pairs:
		rem1 = ((A[1]+1)%n, (A_hat[0]-1+n)%n)
		rem2 = ((A_hat[1]+1)%n, (A[0]-1+n)%n)
		if (rem1 in reflect_square_factor_starts[rem1[0]]) and (rem2 in reflect_square_factor_starts[rem2[0]]):
			return [A, rem1, A_hat, rem2]
	for f in reflect_square_factors:
		f_len = f[1] - f[0] + 1 + n * (f[1] < f[0])
		for of in reflect_square_factor_starts[(f[1]+1)%n]:
			of_len = of[1] - of[0] + 1 + n * (of[1] < of[0])
			if f_len + of_len == n:
				return [f, of]
	return None

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_type_2_half_turn_reflection_tiling(P):
	n = len(P)
	palin_factors = admissible_rotadrome_factors(P, 180)
	reflect_factor_pairs = {}
	for theta in [-45, 0, 45, 90]:
		reflect_factor_pairs[theta] = admissible_gapped_reflect_square_factor_pairs(P, theta)
		# Double up, both for tips and for f1, cf1 iterating since B, refl(B) are *not* symmetric
		reflect_factor_pairs[theta] = (reflect_factor_pairs[theta] + 
			[(fp[1], fp[0]) for fp in reflect_factor_pairs[theta]]) 
	for theta1, theta2 in [(0, 90), (90, 0), (-45, 45), (45, -45)]:
		reflect_factor_tips = {}
		for i in xrange(n):
			for j in range(i) + range(i+1, n):
				reflect_factor_tips[(i, j)] = []
		for f, cf in reflect_factor_pairs[theta2]:
			reflect_factor_tips[(f[0], cf[1])].append((f, cf))			
		# Non-empty B, refl(B)
		for f1, cf1 in reflect_factor_pairs[theta1]:
			dpi1 = ((f1[1]+1)%n, (cf1[0]-1+n)%n)
			dpi2 = ((cf1[1]+1)%n, (f1[0]-1+n)%n)
			# Empty D, refl(D)
			if dpi1 in palin_factors and dpi2 in palin_factors:
				return [f1, dpi1, cf1, dpi2]	
			# Non-empty D, refl(D)
			for f2, cf2 in reflect_factor_tips[(dpi1[0], dpi2[1])]:
				rem1f = None
				if (f2[1]+1)%n == cf1[0]: # |A| = 0 
					rem1f = [f2]
				rem1i = ((f2[1]+1)%n, dpi1[1])
				if rem1f == None and rem1i in palin_factors:
					rem1f = [f2, rem1i]

				rem2f = None
				if cf1[1] == (cf2[0]-1+n)%n: # |C| = 0
					rem2f = [cf2]
				rem2i = (dpi2[0], (cf2[0]-1+n)%n)
				if rem2f == None and rem2i in palin_factors:
					rem2f = [rem2i, cf2]				

				if rem1f and rem2f:
					return [f1] + rem1f + [cf1] + rem2f
		# Empty B, refl(B): D refl(D) A C with A, C, palindromes
		for i in xrange(n):
			for f2, cf2 in reflect_factor_tips[((i+1)%n, i)]:	
				if (f2[1]+1)%n == cf2[0] and (cf2[1]+1)%n == f2[0]:
					return [f2, cf2]
				for p1 in palin_factors:
					if p1[0] != (cf2[1]+1)%n:	 
						continue
					if (p1[1]+1)%n == f2[0]:
						return [f2, cf2, p1]
					for p2 in palin_factors:
						if p2[0] != (p1[1]+1)%n or (p2[1]+1)%n != f2[0]:
							continue
						return [f2, cf2, p1, p2]
	return None

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_type_2_reflection_tiling(P):
	n = len(P)
	mirror_factors = admissible_mirror_factors(P)
	for theta in [-45, 0, 45, 90]:
		reflect_factor_tips = {}
		for i in xrange(n):
			for j in range(i) + range(i+1, n):
				reflect_factor_tips[(i, j)] = []
		for f, cf in admissible_gapped_reflect_square_factor_pairs(P, theta): 
			reflect_factor_tips[(f[0], cf[1])].append((f, cf))
			reflect_factor_tips[(cf[0], f[1])].append((cf, f))
		for A in mirror_factors:
			A_len = A[1] - A[0] + 1 + n * (A[1] < A[0])
			rem1 = ((A[1]+1)%n, (A[0]-1+n)%n)
			for f1, cf1 in reflect_factor_tips[rem1]:
				rem2 = ((f1[1]+1)%n, (cf1[0]-1+n)%n)
				rem2_len = rem2[1] - rem2[0] + 1 + n * (rem2[1] < rem2[0])
				if rem2_len == A_len:
					return [A, f1, rem2, cf1]
				for f2, cf2 in reflect_factor_tips[rem2]:
					rem3 = ((f2[1]+1)%n, (cf2[0]-1+n)%n)
					rem3_len = rem3[1] - rem3[0] + 1 + n * (rem3[1] < rem3[0])
					if rem3_len == A_len:
						return [A, f1, f2, rem3, cf2, cf1]
		for i in xrange(n):
			for f1, cf1 in reflect_factor_tips[((i+1)%n, i)]:	
				rem2 = ((f1[1]+1)%n, (cf1[0]-1+n)%n)
				# Can't be done yet
				for f2, cf2 in reflect_factor_tips[rem2]:
					if (f2[1]+1)%n == cf2[0]:
						return [f1, f2, cf2, cf1]
	return None		

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_type_1_half_turn_reflection_tiling(P):
	n = len(P)
	mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P)
	palin_factors = admissible_rotadrome_factors(P, 180)
	reflect_square_factors = admissible_reflect_square_factors(P)

	palindrome_factor_starts = [[] for i in xrange(n)]
	palindrome_factor_ends = [[] for i in xrange(n)]
	for f in palin_factors:
		palindrome_factor_starts[f[0]].append(f)
		palindrome_factor_ends[f[1]].append(f)
	reflect_square_starts = [[] for i in xrange(n)]
	for f in reflect_square_factors:
		reflect_square_starts[f[0]].append(f)
	
	def is_double_palindrome(F):
		F_len = F[1] - F[0] + 1 + n * (F[1] < F[0])
		for F1 in palindrome_factor_starts[F[0]]:
			F1_len = F1[1] - F1[0] + 1 + n * (F1[1] < F1[0])
			if F1_len == F_len:
				return [F1]
			for F2 in palindrome_factor_ends[F[1]]:
				F2_len = F2[1] - F2[0] + 1 + n * (F2[1] < F2[0])
				if F_len == F1_len + F2_len:
					return [F1, F2]
		return None

	# Non-empty mirrors
	for A, A_hat in mirror_factor_pairs:
		dpi1 = ((A[1]+1)%n, (A_hat[0]-1+n)%n)
		dpi2 = ((A_hat[1]+1)%n, (A[0]-1+n)%n)
		dp1 = is_double_palindrome(dpi1)
		# Other variant (dpi2 is the double palindrome) checked
		# later by symmetry of A, A_hat
		if dp1 and (dpi2 in reflect_square_starts[dpi2[0]]):
			return [A] + dp1 + [A_hat] + [dpi2]

	# Empty mirrors
	for dp1 in reflect_square_factors:
		dp2 = is_double_palindrome(((dp1[1]+1)%n, (dp1[0]-1+n)%n))
		if dp2:
			return [dp1] + dp2
	return None


# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_quarter_turn_tiling(P):
	n = len(P)
	palin_factors = admissible_rotadrome_factors(P, 180)
	ninety_factors = admissible_rotadrome_factors(P, 90)
	ninety_factor_starts = [[] for i in xrange(n)]
	for f in ninety_factors:
		ninety_factor_starts[f[0]].append(f)
	# Factorizations with non-empty palindrome factor
	for C in palin_factors:
		C_len = C[1] - C[0] + 1 + n * (C[1] < C[0])
		for A in ninety_factor_starts[(C[1]+1)%n]:
			A_len = A[1] - A[0] + 1 + n * (A[1] < A[0])
			# One empty 90-drome factor
			# AW: Is this actually possible for polyominoes?
			if A_len + C_len == n: 
				return [C, A]
			# No empty 90-drome factors
			for B in ninety_factor_starts[(A[1]+1)%n]:
				B_len = B[1] - B[0] + 1 + n * (B[1] < B[0])
				if A_len + B_len + C_len == n:
					return [C, A, B]
	# Factorizations with empty palindrome factor
	for A in ninety_factors:
		for B in ninety_factor_starts[(A[1] + 1)%n]:
			AB_len = A[1] - A[0] + 1 + n * (A[1] < A[0]) + B[1] - B[0] + 1 + n * (B[1] < B[0])
			if AB_len == n:
				return [A, B]
	return None

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_translation_tiling(P):
	n = len(P)
	factors = admissible_mirror_factors(P)
	factor_starts = [[] for i in xrange(n)]
	factor_ends = [[] for i in xrange(n)]
	for f in factors:
		factor_starts[f[0]].append(f)
		factor_ends[f[1]].append(f)
	for i in xrange(n):
		for A in factor_starts[i]:
			for B in factor_starts[(A[1]+1)%n]:
				AB_len = A[1] - A[0] + 1 + n * (A[1] < A[0]) + B[1] - B[0] + 1 + n * (B[1] < B[0])
				if AB_len > n/2:
					continue
				if AB_len == n/2:
					return [A, B]
				C = ((B[1]+1)%n, (A[0]+n/2-1+n)%n) 
				if C in factor_starts[C[0]]:
					return [A, B, C]
	return None

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_half_turn_tiling(P):
	n = len(P)
	mirror_factor_pairs = admissible_gapped_mirror_factor_pairs(P)
	palindrome_factor_starts = [[] for i in xrange(n)]
	palindrome_factor_ends = [[] for i in xrange(n)]
	for f in admissible_rotadrome_factors(P, 180):
		palindrome_factor_starts[f[0]].append(f)
		palindrome_factor_ends[f[1]].append(f)

	def is_double_palindrome(F):
		F_len = F[1] - F[0] + 1 + n * (F[1] < F[0])
		for F1 in palindrome_factor_starts[F[0]]:
			F1_len = F1[1] - F1[0] + 1 + n * (F1[1] < F1[0])
			if F1_len == F_len:
				return [F1]
			for F2 in palindrome_factor_ends[F[1]]:
				F2_len = F2[1] - F2[0] + 1 + n * (F2[1] < F2[0])
				if F_len == F1_len + F2_len:
					return [F1, F2]
		return None

	for A, A_hat in mirror_factor_pairs:
		dpi1 = ((A[1]+1)%n, (A_hat[0]-1+n)%n)
		dpi2 = ((A_hat[1]+1)%n, (A[0]-1+n)%n)
		dp1 = is_double_palindrome(dpi1)
		dp2 = is_double_palindrome(dpi2)
		if dp1 and dp2:
			return [A] + dp1 + [A_hat] + dp2

	for i in xrange(n):
		for j in range(i) + range(i+1, n):
			dpi1 = (i, (j-1+n)%n)
			dpi2 = (j%n, (i-1+n)%n)
			dp1 = is_double_palindrome(dpi1)
			dp2 = is_double_palindrome(dpi2)
			if dp1 and dp2:
				return [A] + dp1 + [A_hat] + dp2
	return None

# Input: a polyomino boundary word
# Output: a witness boundary decomposition or None
def has_isohedral_tiling(P):
	return (has_half_turn_tiling(P) or has_translation_tiling(P) or has_quarter_turn_tiling(P) or 
		has_type_1_reflection_tiling(P) or has_type_2_reflection_tiling(P) or
		has_type_1_half_turn_reflection_tiling(P) or has_type_2_half_turn_reflection_tiling(P))

class TestStuff(unittest.TestCase):

	def setUp(self):
		pass
	
	def test__admissible_rotadrome_factors(self):
		self.assertEqual(set(admissible_rotadrome_factors(['N', 'E', 'S', 'W'], 180)), 
			set([(0, 0), (1, 1), (2, 2), (3, 3)])) 

	def test__admissible_reflect_square_factors(self):
		self.assertEqual(set(admissible_reflect_square_factors(['N', 'E', 'E', 'S', 'W', 'W'])),
			set([(0, 1), (2, 3), (3, 4), (5, 0), (1, 2), (4, 5)]))
		self.assertEqual(set(admissible_reflect_square_factors(['N', 'N', 'E', 'E', 'S', 'S', 'W', 'W'])),
			set([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0), 
				(0, 3), (2, 5), (4, 7), (6, 1)])) 

	def test__admissible_gapped_reflect_square_factor_pairs(self):
		E = [(i, i) for i in xrange(4)]
		self.assertEqual(set(admissible_gapped_reflect_square_factor_pairs(['N', 'E', 'S', 'W'], -45)),
			set([(E[0], E[3]), (E[1], E[2])]))
		self.assertEqual(set(admissible_gapped_reflect_square_factor_pairs(['N', 'E', 'S', 'W'], 45)),
			set([(E[0], E[1]), (E[2], E[3])])) 
		self.assertEqual(set(admissible_gapped_reflect_square_factor_pairs(['N', 'E', 'S', 'W'], 0)),
			set([(E[0], E[2])])) 
		self.assertEqual(set(admissible_gapped_reflect_square_factor_pairs(['N', 'E', 'S', 'W'], 90)),
			set([(E[1], E[3])])) 

	def test__has_type_1_reflection_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_type_1_reflection_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 2xi rectangles
		for i in [1,3,4,5]:
			self.assertTrue(has_type_1_reflection_tiling(['N', 'N'] + ['E']*i + ['S', 'S'] + ['W']*i))
		# 3xi rectangles
		for i in [4,5,6,7]:
			self.assertTrue(has_type_1_reflection_tiling(['N', 'N'] + ['E']*i + ['S', 'S'] + ['W']*i))
		# Small tetris pieces
		self.assertFalse(has_type_1_reflection_tiling(['N', 'W', 'N', 'E', 'N', 'E', 'N', 'E', 'S', 'E',
			'S', 'E', 'S', 'W', 'S', 'W', 'N', 'W', 'S', 'W']))	
		self.assertTrue(has_type_1_reflection_tiling([
			'N', 'N', 'W', 'N', 'E', 'N', 
			'N', 'N', 'E', 'N', 'W', 'N', 
			'E', 'E', 'E', 'E', 'N', 'E', 'S', 'E', 'E', 'E', 
			'S', 'E', 'S', 'E', 'S', 'E',
			'S', 'W', 'S', 'W', 'S', 'W', 
			'W', 'W', 'W', 'N', 'W', 'S', 'W', 'W', 'W', 'W']))

	def test__has_type_2_reflection_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_type_2_reflection_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 2xi rectangles
		for i in [1,3,4,5]:
			self.assertTrue(has_type_2_reflection_tiling(['N', 'N'] + ['E']*i + ['S', 'S'] + ['W']*i))
		# 3xi rectangles
		for i in [4,5,6,7]:
			self.assertTrue(has_type_2_reflection_tiling(['N', 'N'] + ['E']*i + ['S', 'S'] + ['W']*i))
		# Tetris pieces (non-empty A, A hat)
		self.assertTrue(has_type_2_reflection_tiling([
			'N', 'N', 'E', 'N', 'W', 'N', 	
			'N', 'N', 'E', 'E', 'N', 'N',
			'E', 'E', 'E', 'E', 'N', 'E', 'E', 'S', 'E',
			'S', 'S', 'W', 'W', 'S', 'S',
			'S', 'E', 'S', 'W', 'S', 'S',
			'W', 'N', 'W', 'W', 'S', 'W', 'W', 'W', 'W']))
		self.assertFalse(has_type_2_reflection_tiling([
			'N', 'N', 'E', 'N', 'W', 'N', 	
			'N', 'N', 'E', 'E', 'N', 'N',
			'E', 'E', 'E', 'E', 'N', 'E', 'E', 'S', 'E',
			'S', 'S', 'W', 'W', 'S', 'S',
			'S', 'W', 'S', 'E', 'S', 'S',
			'W', 'N', 'W', 'W', 'S', 'W', 'W', 'W', 'W']))
		# Tetris pieces (empty A, A hat)
		self.assertTrue(has_type_2_reflection_tiling([
			'W', 'N', 'W', 'N', 
			'N', 'N', 'E', 'E',
			'S', 'S', 'E', 'E',
			'W', 'S', 'W', 'S'])) 
		self.assertFalse(has_type_2_reflection_tiling([
			'W', 'N', 'W', 'N', 
			'N', 'N', 'E', 'E',
			'S', 'S', 'E', 'E',
			'W', 'W', 'S', 'S'])) 

	def test__has_type_2_half_turn_reflection_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_type_2_half_turn_reflection_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# Tetris pieces (non-empty D, non-empty B, refl(B), empty A, C)
		self.assertTrue(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'E', 'S', 'W', 'S',
			'W', 'N', 'W', 'S', 'W']))
		self.assertFalse(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'W', 'S', 'E', 'S',
			'W', 'N', 'W', 'S', 'W']))
		# Tetris pieces (non-empty D, non-empty B, refl(B), non-empty A, C)
		self.assertTrue(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'W', 'N', 'E', 'N',
			'N', 'E', 'N', 'W', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'E', 'S', 'W', 'S',
			'S', 'S', 'S', 'S',
			'W', 'N', 'W', 'S', 'W']))
		self.assertFalse(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'W', 'N', 'E', 'N',
			'N', 'E', 'N', 'W', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'E', 'S', 'W', 'S',
			'S', 'E', 'S', 'W', 'S', 'S',
			'W', 'N', 'W', 'S', 'W']))
		# Tetris pieces (empty D, non-empty B, refl(B), non-empty A, C)
		self.assertTrue(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'W', 'N', 'E', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'S', 'S', 'S',
			'W', 'N', 'W', 'S', 'W']))
		self.assertFalse(has_type_2_half_turn_reflection_tiling([
			'N', 'E', 'N', 'W', 'W', 'N', 'E', 'N',
			'E', 'N', 'E', 'S', 'E',
			'S', 'S', 'E', 'S', 'W', 'S',
			'W', 'N', 'W', 'S', 'W']))
		

	def test__has_type_1_half_turn_reflection_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_type_1_half_turn_reflection_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# Tetris pieces (non-empty A, A hat)
		self.assertFalse(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'W', 'N', 'E', 'N', 
			'N', 'N', 'E', 'N', 'W', 'N', 
			'E', 'E', 'E', 'E', 'N', 'E', 'S', 'E', 'E', 'E', 
			'S', 'S', 'S', 'S', 
			'S', 'E', 'S', 'W', 'S', 'S', 
			'W', 'W', 'W', 'N', 'W', 'S', 'W', 'W', 'W', 'W']))
		self.assertTrue(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'N', 'W', 'N', 'E', 'N', 
			'N', 'N', 'N', 'E', 'N', 'W', 'N', 
			'E', 'E', 'E', 'E', 'N', 'E', 'S', 'E', 'E', 'E', 
			'S', 'S', 'S', 'S', 
			'S', 'E', 'S', 'W', 'W', 'S', 'E', 'S',	
			'W', 'W', 'W', 'N', 'W', 'S', 'W', 'W', 'W', 'W']))
		# Tetris pieces (empty A, A hat)
		self.assertTrue(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'N', 'N', 'N', 'N',
			'E', 'S', 'E', 'S', 'E', 'S',
			'W', 'S', 'W', 'W', 'S', 'W']))
		self.assertTrue(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'W', 'N', 'E', 'E', 'N', 'W', 'N', 'N',
			'E', 'S', 'E', 'S', 'E', 'S',
			'W', 'S', 'W', 'W', 'S', 'W']))
		self.assertFalse(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'W', 'N', 'E', 'N', 'N', 'N',
			'E', 'S', 'E', 'S', 'E', 'S',
			'W', 'S', 'W', 'W', 'S', 'W']))
		self.assertFalse(has_type_1_half_turn_reflection_tiling([
			'N', 'N', 'W', 'N', 'E', 'E', 'N', 'W', 'N', 'N',
			'E', 'S', 'S', 'E', 'E', 'S',
			'W', 'S', 'W', 'W', 'S', 'W']))

	def test__has_quarter_turn_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_quarter_turn_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 1xi rectangles
		self.assertTrue(has_quarter_turn_tiling(['N'] + ['E']*2 + ['S'] + ['W']*2))
		for i in [3,4,5]:
			self.assertFalse(has_quarter_turn_tiling(['N'] + ['E']*i + ['S'] + ['W']*i))
		# Small tetris pieces
		self.assertTrue(has_quarter_turn_tiling(['N', 'N', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'W']))
		self.assertFalse(has_quarter_turn_tiling(['N', 'N', 'E', 'N', 'E', 'S', 'E', 'S', 'S', 'W', 'N', 'W', 'S', 'W']))

	def test__has_translation_tiling(self):
		# Squares of various sizes
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N']*i + ['E']*i + ['S']*i + ['W']*i))
		# 1xi, ix1 rectangles
		for i in [1,2,3,4,5]:
			self.assertTrue(has_translation_tiling(['N'] + ['E']*i + ['S'] + ['W']*i))
			self.assertTrue(has_translation_tiling(['N']*i + ['E'] + ['S']*i + ['W']))
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
			self.assertTrue(has_half_turn_tiling(['N']*i + ['E'] + ['S']*i + ['W']))
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
	exit(0)
		


