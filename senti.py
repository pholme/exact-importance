# this code was written by Petter Holme spring / summer 2017
# it calculates the time to discovery or extinction in the SIR model given a set of sentinels
# the outbreak size is averaged over all nodes as seed nodes
# for more detailed comments, see the infmax.py program

import networkx as nx
from sys import argv
from sympy.abc import x
from sympy import Poly
from gc import collect

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# returns a set of nodes that could be infected and one set that could recover

def get_infect_reco ():
	global G

	infect = {}
	for u,v in G.edges_iter():
		if G.node[u]['state'] == 'I' and (G.node[v]['state'] == 'S' or G.node[v]['state'] == 'SENTI'):
			if v in infect:
				infect[v] += 1
			else:
				infect[v] = 1
		elif (G.node[u]['state'] == 'S' or G.node[u]['state'] == 'SENTI') and G.node[v]['state'] == 'I':
			if u in infect:
				infect[u] += 1
			else:
				infect[u] = 1

	return infect, [v for v in G.nodes() if G.node[v]['state'] == 'I']

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def stepdown (wnum0,wden0,dtnum0,dtden0):
	global G, tnum, tden

	wnum = wnum0.mul_ground(1)
	wden = wden0.mul_ground(1)
	dtnum = dtnum0.mul_ground(1)
	dtden = dtden0.mul_ground(1)

	infectables, recoverables = get_infect_reco()
	si = sum(infectables.values())
	sr = len(recoverables)

	if si == 0 and sr == 0:
		a = wnum.mul(dtnum)
		b = wden.mul(dtden)
		c = b.mul(tden)
		a = a.mul(tden)
		b = b.mul(tnum)
		a = a.add(b)
		b = a.gcd(c)
		tnum,no = a.div(b)
		tden,no = c.div(b)
		
		collect()
		return
	
	den = Poly(si * x + sr, x)
	a = den.mul(dtnum)
	a = a.add(dtden)
	b = den.mul(dtden)
	c = a.gcd(b)
	dtnum,no = a.div(c)
	dtden,no = b.div(c)
	
	nden = wden.mul(den)
	for you, num in sorted(infectables.iteritems()):
		nnum = wnum.mul(Poly(num * x,x))

		if G.node[you]['state'] == 'SENTI':
			a = nnum.mul(dtnum)
			b = nden.mul(dtden)
			c = b.mul(tden)
			a = a.mul(tden)
			b = b.mul(tnum)
			a = a.add(b)
			b = a.gcd(c)
			tnum,no = a.div(b)
			tden,no = c.div(b)
		else:
			G.node[you]['state'] = 'I'
			stepdown(nnum,nden,dtnum,dtden)
			G.node[you]['state'] = 'S'
	for you in sorted(recoverables):
		G.node[you]['state'] = 'R'	
		stepdown(wnum,nden,dtnum,dtden)
		G.node[you]['state'] = 'I'
	collect()

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

if __name__ == "__main__":
	global G, tnum, tden

	if len(argv) < 2:
		print 'usage python senti.py [# links] [links] <sentinels>'
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
		G0.node[me]['state'] = 'SENTI'

	G = nx.convert_node_labels_to_integers(G0,label_attribute='id')

	stnum = Poly(0,x)
	stden = Poly(1,x)
	for v in G.nodes():
		for u in G.nodes():
			if G.node[u]['state'] != 'SENTI':
				G.node[u]['state'] = 'S'
		if G.node[v]['state'] != 'SENTI':
			G.node[v]['state'] = 'I'
		tnum = Poly(0,x)
		tden = Poly(1,x)

		stepdown(Poly(1,x),Poly(1,x),Poly(0,x),Poly(1,x))
		
		a = tnum.mul(stden)
		b = stnum.mul(tden)
		a = a.add(b)
		b = tden.mul(stden)
		c = a.gcd(b)
		stnum,no = a.div(c)
		stden,no = b.div(c)
	
	a = stden.mul_ground(G.number_of_nodes())
	b = stnum.gcd(a)
	c,no = stnum.div(b)
	d,no = a.div(b)
	print s.strip() + ', (' + str(c.as_expr()) + ')/(' + str(d.as_expr()) + ')'

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
