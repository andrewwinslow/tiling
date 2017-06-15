
import polyomino 

# Input: an input filename, grid dimensions, and an output filename 
# Output: none (a postscript file named outfile is created and filled with a
# grid of drawings of the polyomino boundary words contained in the file named 
# infile)
def draw_polyominoes(infile, rows, cols, outfile): 

	def polyominotops(path, orig = (300,300), scale = 30): 
		ox,oy = orig
		s = "%f %f newpath moveto\n"%(ox,oy)
		for d in list(path):
			ox += polyomino.dir2vec[d][0]*scale
			oy += polyomino.dir2vec[d][1]*scale
			s += "%f %f lineto\n"%(ox,oy)
		s += "stroke\n"
		return s

	scale = 2
	f = open(infile, 'r')
	outfile = open(outfile, 'w')
	outfile.write("%!PS\n\n")
	i = 0
	j = 0
	for line in f:
		orig = (50+((i % cols)*19*scale, 30+(j % rows)*19*scale)
		outfile.write(polyominotops(line[:-1], orig, scale)) 
		i += 1
		if (i == cols)
			j += 1
		if (j == rows)
			outfile.write("showpage\n")
	infile.close()
	outfile.close()
	
draw_polyominoes("dicube_unfolding_boundary_words.txt", 10, 10, "tilings.ps")

