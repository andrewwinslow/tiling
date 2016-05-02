
import polyomino
import isohedral
import polycube
import sys
import random

uf_limit = 1000
ssuf_limit = 1000

polycube_count = 0
for P in polycube.enumerate_polycubes(7):
	polycube_count = polycube_count + 1
	print ""
	print "Processing",
	for c in P:
		print str(c),
	print "(heptacube #" + str(polycube_count) + ")..."
	success = False
	seed = 20
	while not success:
		random.seed(seed)
		uf_count = 0
		ssuf_count = 0
		for uf in polycube.unfoldings(P):
			if uf_count > uf_limit or ssuf_count > ssuf_limit:
				break
			uf_count = uf_count + 1
			cuf = polyomino.cancel(uf)
			if not polyomino.is_polyomino(cuf):
				continue
			ssuf_count = ssuf_count + 1
			soln = isohedral.has_half_turn_tiling(cuf)
			if soln:
				print ""
				print uf
				print cuf
				print soln 
				success = True
				break	
		if not success:
			print ""
		seed = seed + 1


