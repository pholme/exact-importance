// this file contains get_inf_rec making a list of vertices that can be
// infected or recovered, using the IGRAPH implementation of the VF2 graph
// isomorphims algo

#include "poly.h"
#include <igraph.h>

extern GLOBAL g;
extern NODE n[MAXN];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// This routine makes an igraph representation of the network disregarding 
// R nodes

void make_igraph (igraph_t *a, igraph_vector_int_t *c) {
	int i, j, u, v;

	// First do a simple BFS to eliminate isolated components without I
	// (they will not contribute to the result). This part doesn't speed
	// things more than 1-2% and could maybe be eliminated for clarity.

	for (i = g.nn = 0; i < g.n; i++) if (IS_INFECTIOUS(i)) {
		n[i].state += 10; // set a temporary state for nodes in the igraph
		// which are also nodes marked as checked in the BFS
		g.check[g.nn++] = i;
	}

	for (i = 0; i < g.nn; i++) {
		u = g.check[i];
		for (j = 0; j < n[u].deg; j++) {
			v = n[u].nb[j];
			if (IS_SUSCEPTIBLE(v)) {
				n[v].state += 10;
				g.check[g.nn++] = v;
			}
		}
	}

	for (i = j = 0; i < g.n; i++) {
		if (n[i].state >= 10) {
			n[i].id = j++;
			n[i].state -= 10; // reset the state
		} else n[i].id = -1; // mark outsiders as negative
	}

	// then create and fill up the graph

	igraph_empty(a, g.nn, IGRAPH_UNDIRECTED);
	igraph_vector_int_init(c, g.nn);

	for (i = 0; i < g.n; i++) {
		u = n[i].id;
		if (u >= 0) {
			VECTOR(*c)[u] = n[i].state;
			for (j = 0; j < n[i].deg; j++) {
				v = n[i].nb[j];
				if (i < v) if (n[v].id >= 0) igraph_add_edge(a, u, n[v].id);
			}
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// is it equivalent if me or you becomes infected?

int iso_i (int me, int you) {
	igraph_bool_t isomorphic;
	igraph_t gme;
	igraph_vector_int_t cme, cyou;

	// create an igraph
	make_igraph(&gme, &cme);

	igraph_vector_int_init(&cyou, igraph_vector_int_size(&cme));
	igraph_vector_int_update(&cyou, &cme);

	VECTOR(cme)[n[me].id] = VECTOR(cyou)[n[you].id] = INFECTIOUS;

	igraph_isomorphic_vf2(&gme, &gme, &cme, &cyou, NULL, NULL,
			&isomorphic, NULL, NULL, NULL, NULL, NULL);

	igraph_destroy(&gme);
	igraph_vector_int_destroy(&cme);
	igraph_vector_int_destroy(&cyou);

	if (isomorphic) g.niso_i++;
	else g.nniso_i++;

	return (int) isomorphic;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// is it equivalent if me or you recovers?

int iso_r (int me, int you) {
	igraph_bool_t isomorphic;
	igraph_t gme, gyou;
	igraph_vector_int_t cme, cyou;

	// create two igraphs corresponding to me and you being infected
	n[me].state = RECOVERED;
	make_igraph(&gme, &cme);
	n[me].state = INFECTIOUS;

	n[you].state = RECOVERED;
	make_igraph(&gyou, &cyou);
	n[you].state = INFECTIOUS;

	if (igraph_vcount(&gme) == igraph_vcount(&gyou)) 
		igraph_isomorphic_vf2(&gme, &gyou, &cme, &cyou, NULL, NULL,
			&isomorphic, NULL, NULL, NULL, NULL, NULL);
	else isomorphic = FALSE;

	igraph_destroy(&gme);
	igraph_destroy(&gyou);
	igraph_vector_int_destroy(&cme);
	igraph_vector_int_destroy(&cyou);

	if (isomorphic) g.niso_r++;
	else g.nniso_r++;

	return (int) isomorphic;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// making a list of nodes that can be infected or recovered

void get_inf_rec (int *ninfectables, int *to_step_to, int *nto_step_to,
		int *multiplicity, int *nrecoverables, int *nsi) {
	int i, j, me, you;

	*ninfectables = *nto_step_to = *nrecoverables = *nsi = 0;

	for (i = 0; i < g.n; i++) multiplicity[i] = 0;

	for (me = 0; me < g.n; me++) {
		if (IS_SUSCEPTIBLE(me)) {
			for (i = 0; i < n[me].deg; i++) {
				you = n[me].nb[i];
				if (IS_INFECTIOUS(you)) {
					multiplicity[me]++;
					(*nsi)++;
				}
			}
			if (multiplicity[me] > 0) to_step_to[(*ninfectables)++] = me;
		} else {
			if (IS_INFECTIOUS(me)) {
				multiplicity[me] = 1;
				g.recoverables[(*nrecoverables)++] = me;
			}
		}
	}

	// eliminate nodes that whose infection leads to equivalent states
	for (i = 0; i < (*ninfectables) - 1; i++) {
		me = to_step_to[i];
		if (multiplicity[me] > 0) {
			for (j = i + 1; j < (*ninfectables); j++) {
				you = to_step_to[j];
				if (iso_i(me, you)) { // merge if equivalent
					multiplicity[me] += multiplicity[you];
					multiplicity[you] = 0;
				}
			}
		}
	}

	for (i = *nto_step_to = 0; i < (*ninfectables); i++) {
		me = to_step_to[i];
		if (multiplicity[me] > 0) to_step_to[(*nto_step_to)++] = me;
	}
	*ninfectables = *nto_step_to;

	// eliminate nodes that whose recovery leads to equivalent states
	for (i = 0; i < (*nrecoverables) - 1; i++) {
		me = g.recoverables[i];
		if (multiplicity[me] > 0) {
			for (j = i + 1; j < (*nrecoverables); j++) {
				you = g.recoverables[j];
				if (iso_r(me, you)) { // merge if equivalent
					multiplicity[me] += multiplicity[you];
					multiplicity[you] = 0;
				}
			}
		}
	}

	for (i = 0; i < (*nrecoverables); i++) {
		me = g.recoverables[i];
		if (multiplicity[me] > 0) to_step_to[(*nto_step_to)++] = me;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
