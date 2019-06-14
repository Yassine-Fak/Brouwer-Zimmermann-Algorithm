# A simple toy ISD
# PK - 2018


def find_lo_weight(C,maxiter=1000000,wtarget=1,wcomb=2,nthreads=1,niterperthread=1):
	from sage.combinat.gray_codes import combinations
	import time

	if GF(2) != C.base_field():
		raise NotImplementedError

	n,k = C.length(), C.dimension()
	G = C.generator_matrix()
	supp = range(n)
	minwd = G.row(0) # arbitrary
	minw = minwd.hamming_weight()
	it = 0
	nthreads = min(maxiter/niterperthread, nthreads)
	epoch = time.clock()

	@parallel(nthreads)
	def do_Iset(curminw,nbiter):
		set_random_seed()
		curminwd = None
		for i in xrange(nbiter):
			while True:
				Iset = sample(supp,k)
				Gis = G.matrix_from_columns(Iset)
				try:
					Gis_inv = Gis.inverse()
					break
				except ZeroDivisionError:
					continue
			Glw = Gis_inv * G
			# weight one is simple
			for c in xrange(k):
				cc = Glw.row(c)
				cw = cc.hamming_weight()
				if cw < curminw:
					curminwd,curminw = cc,cw
			# further weights wiz Gray codes
			for w in range(wcomb)[1:]:
				base = [1]*w
				base.extend([0]*(k-w))
				v = vector(base)
				cc = v*Glw
				cw = cc.hamming_weight() #ugly repeat 1
				if cw < curminw:
					curminwd,curminw = cc,cw
				for i,j in combinations(k,w):
					cc += Glw.row(i)
					cc += Glw.row(j)
					cw = cc.hamming_weight()
					if cw < curminw: #ugly repeat 2
						curminwd,curminw = cc,cw
		return (curminw,curminwd)

	while it < maxiter:
		if 1 == nthreads:
			(cw,cc) = do_Iset(minw,niterperthread)
		else:
			run = do_Iset([(minw,niterperthread)]*nthreads)
			res = list(run)
			(cw,cc) = min(res)[1]
		if cw < minw:
			minwd,minw = cc,cw
			print "Found a new codeword of weight "+str(cw)+" @ time "+str(time.clock() - epoch)
		# done once per iteration, as not expected to be successful many times?
		if minw <= wtarget:
			return minwd
		it += nthreads*niterperthread
	return minwd

# TODO always include the error in the linear combs
def find_lo_error(C,w,maxiter=1000000,maxerrweight=1,wcomb=2,nthreads=1,niterperthread=1):
	s = C.syndrome(matrix(w).transpose())
	if (0 == vector(s).hamming_weight()):
		return (w,w+w)

	G = C.generator_matrix()
	Gg = block_matrix([G,matrix(w)],nrows=2)
	Cc = LinearCode(Gg)

	e = find_lo_weight(Cc,maxiter,maxerrweight,wcomb,nthreads,niterperthread)
	wc = w+e
	s = C.syndrome(matrix(wc).transpose())
	if (0 == vector(s).hamming_weight()):
		print "Found a codeword with distance "+str(e.hamming_weight())
		return (wc,e)
	else:
		print "Failure"
		return (None,None)
