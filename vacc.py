# this code was written by Petter Holme spring / summer 2017
# it calculates the outbreak size of the SIR model given a set of nodes deleted
# the outbreak size is averaged over all nodes as seed nodes
# for more detailed comments, see the infmax.py program

import networkx as nx
from sys import argv
from sympy.abc import x
from sympy import Poly

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def get_infect_reco ():
	global G

	infect = {}
	for u,v in G.edges_iter():
		if G.node[u]['state'] == 'I' and G.node[v]['state'] == 'S':
			if v in infect:
				infect[v] += 1
			else:
				infect[v] = 1
		elif G.node[u]['state'] == 'S' and G.node[v]['state'] == 'I':
			if u in infect:
				infect[u] += 1
			else:
				infect[u] = 1

	return infect, [v for v in G.nodes() if G.node[v]['state'] == 'I']

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def obsize (G):
		
	return len([v for v in G.nodes() if G.node[v]['state'] == 'R' or G.node[v]['state'] == 'I'])

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def stepdown (wnum0,wden0):
	global G, onum, oden

	wnum = wnum0.mul_ground(1)
	wden = wden0.mul_ground(1)

	infectables, recoverables = get_infect_reco()
	si = sum(infectables.values())

	if si == 0:
		a = onum.mul(wden)
		b = wnum.mul(oden).mul_ground(obsize(G))
		c = oden.mul(wden)
		d = a.add(b)
		a = d.gcd(c)
		onum,no = d.div(a)
		oden,no = c.div(a)
		return

	sr = len(recoverables)
	den = Poly(si * x + sr, x)
	
	for you, num in infectables.iteritems():
		G.node[you]['state'] = 'I'
		
		stepdown(wnum.mul(Poly(num * x,x)),wden.mul(den))
		G.node[you]['state'] = 'S'
	for you in recoverables:
		G.node[you]['state'] = 'R'
		
		stepdown(wnum,wden.mul(den))
		G.node[you]['state'] = 'I'

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

if __name__ == "__main__":
	global G, onum, oden

	if len(argv) < 2:
		print 'usage python vacc.py [# links] [links] <vaccinees>'
		exit()

	G0 = nx.Graph()

	nl = int(argv[1])

	for i in range(nl):
		G0.add_edge(argv[2 + 2 * i], argv[3 + 2 * i])

	s = ''
	for v in G0.nodes():
		G0.node[v]['state'] = 'S'
	for i in range(2 + 2 * nl,len(argv)):
		me = argv[i]
		s += me + ' '
		if me not in G0.nodes():
			print me, 'not in G'
			exit()
		G0.node[me]['state'] = 'V'

	G = nx.convert_node_labels_to_integers(G0,label_attribute='id')
	
	sonum = Poly(0,x)
	soden = Poly(1,x)
	for v in G.nodes():
		if G.node[v]['state'] != 'V':
			for u in G.nodes():
				if G.node[u]['state'] != 'V':
					G.node[u]['state'] = 'S'
			G.node[v]['state'] = 'I'
			onum = Poly(0,x)
			oden = Poly(1,x)
			
			stepdown(Poly(1,x),Poly(1,x))

			a = onum.mul(soden)
			b = sonum.mul(oden)
			c = a.add(b)
			d = oden.mul(soden)
			a = d.gcd(c)
			sonum,no = c.div(a)
			soden,no = d.div(a)

			print '(' + str(sonum.as_expr()) + ')/(' + str(soden.as_expr()) + ')'
			print '-------------'
	
	a = soden.mul_ground(G.number_of_nodes())
	b = sonum.gcd(a)
	c,no = sonum.div(b)
	d,no = a.div(b)
	print s.strip() + ', (' + str(c.as_expr()) + ')/(' + str(d.as_expr()) + ')'

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
