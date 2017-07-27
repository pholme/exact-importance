'''

This program analyzes the output of infmax, vacc and senti.
It assumes there are files with the name 'fname' in the directories
res/{vacc,infmax,senti}/ that contains output from all three programs.
These files should contain output from all sets of active nodes. For the
example of one link and n = 1, this would be achieved by.

./infmax 1 0 1 0 >> res/infmax/file
./infmax 1 0 1 1 >> res/infmax/file
./vacc 1 0 1 0 >> res/vacc/file
./vacc 1 0 1 1 >> res/vacc/file
./senti 1 0 1 0 >> res/senti/file
./senti 1 0 1 1 >> res/senti/file

python eval_res.py file

the output shows the intervals of beta with distinct optimal sets of vertices
(both by exact algebraic expressions and numerical values, separated by "=")

'''


from sys import argv
from sympy import S, oo, Interval, Poly, N
from sympy.abc import x
from sympy.polys.polytools import poly_from_expr
from sympy.solvers.inequalities import solve_poly_inequality
from itertools import combinations

tol = 1e-5
senses = ['infmax','vacc','senti']
large = 10000
larger = large + 10000

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def eliminate (s, t):
	global maxsol, maximize

	a = maxsol[s]['num'].mul(maxsol[t]['den'])
	b = maxsol[t]['num'].mul(maxsol[s]['den'])
	p = a.sub(b)

	if maximize:
		a = solve_poly_inequality(p, '>')
	else:
		a = solve_poly_inequality(p, '<')

	# make a union of intervals

	c = []
	for b in a:
		if b.measure > tol:
			if b.sup > 0:
				c.append(b)

	if len(c) == 0:
		maxsol[s]['interval'] = S.EmptySet
		return

	b = c[0]
	for i in range(1,len(c)):
		b = b.union(c[i])

	# delete the suboptimal regions

	a = b.complement(S.Reals)
	maxsol[s]['interval'] = maxsol[s]['interval'].intersection(b)
	maxsol[t]['interval'] = maxsol[t]['interval'].intersection(a)

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# being a bit more generous than just poly_from_expr

def make_poly (s):

	try:
		a,no = poly_from_expr(s)
	except:
		if not 'x' in s:
			s = s.strip('()')
			if s.isdigit():
				return Poly(int(s),x)
		print "couldn't convert", s, 'to polynomial'
	return a

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def different (p, s, ii):

	for m in senses:
		if len(p[m]) == 0:
			return True
		if ii[m] < len(s[m]):
			if set(p[m]) != set(s[m][ii[m]]['nodes']):
				return True
	return False

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def get_sol (d,fname):
	global maximize, maxsol

	f = open('res/' + d + '/' + fname)
	if d == 'infmax':
		maximize = True
	else:
		maximize = False

	polys = {}
	maxsol = {}

	for l in f:
		a = l.split(',')
		s = a[0].strip()
		y = a[1].strip()

		if s not in polys.values(): # ignore duplicates
			if y in polys:
				polys[y].append(s)
			else:
				polys[y] = [s]

				a = y.split('/')
				
				maxsol[y] = {}
				maxsol[y]['num'] = make_poly(a[0])
				if len(a) > 1:
					maxsol[y]['den'] = make_poly(a[1])
				else:
					maxsol[y]['den'] = Poly(1,x)
				maxsol[y]['interval'] = Interval.Ropen(0,oo)

	f.close()

	for s,t in combinations(maxsol,2):
		eliminate(s,t)

	solint = []
	for s in maxsol:
		if maxsol[s]['interval'].is_Union:
			a = list(maxsol[s]['interval'].args)
		else:
			a = [maxsol[s]['interval']]
		for b in a:
			if b.measure > tol:
				solint.append({})
				solint[-1]['nodes'] = polys[s]
				solint[-1]['beginning'] = b.inf

	solint.sort(key = lambda c: c['beginning'])
	
	return solint

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

if __name__ == "__main__":

	if len(argv) < 2:
		print 'usage python eval_res.py [file name]'
		exit(1)

	if '/' in argv[1]:
		print 'just give name of file inside directory'
		exit()
	
	s = {}
	i = {}
	prn = {}
	for m in senses:
		s[m] = get_sol(m,argv[1]) 
		prn[m] = []
		i[m] = 0

	smallest = 0
	while True:
		if different(prn,s,i):
			print str(smallest), '=', str(N(smallest))
			print 'I', str(s['infmax'][i['infmax']]['nodes']) + ', V', str(s['vacc'][i['vacc']]['nodes']) + ', S', str(s['senti'][i['senti']]['nodes'])
			for m in senses:
				prn[m] = s[m][i[m]]['nodes']

		# step up smallest
		smallest = larger
		adv = []
		for m in senses:
			k = i[m] + 1
			if k < len(s[m]):
				a = s[m][k]['beginning']
				if a == smallest:
					adv.append(m)
				elif a < smallest:
					smallest = a
					adv = [m]

		if smallest == larger:
			break

		for m in adv:
			i[m] += 1

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #







	