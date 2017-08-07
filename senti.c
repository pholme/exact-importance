// C code corresponding to senti.py (using the FLINT library)
// see senti.py for more detailed comments

#include "poly.h"

NODE n[N];
GLOBAL g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// making a list of nodes that can be infected or recovered (very naive)

void get_infect_reco (int *infectables, int *ninfectables,
		int *recoverables, int *nrecoverables, int *infways) {
	int i, me, you;

	for (i = 0; i < g.n; i++) infways[i] = 0;

	for (me = 0; me < g.n; me++) {
		if ((n[me].state == S) || (n[me].state == SENTI)) {
			for (i = 0; i < n[me].deg; i++) {
				you = n[me].nb[i];
				if (n[you].state == I) infways[me]++;
			}
			if (infways[me] > 0) infectables[(*ninfectables)++] = me;
		} else {
			if (n[me].state == I) recoverables[(*nrecoverables)++] = me;
		}
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void stepdown (fmpz_poly_t wnum, fmpz_poly_t wden0, fmpz_poly_t donum0, fmpz_poly_t doden0) {
	int i, si = 0, you, infways[N];
	int infectables[N], ninfectables = 0, recoverables[N], nrecoverables = 0;
	fmpz_poly_t wden, donum, doden, nnum;

	fmpz_poly_init(donum);
	fmpz_poly_init(doden);
	fmpz_poly_init(wden);
	fmpz_poly_init(nnum);

	fmpz_poly_set(wden, wden0);
	fmpz_poly_set(donum, donum0);
	fmpz_poly_set(doden, doden0);

	get_infect_reco(infectables, &ninfectables, recoverables, &nrecoverables, infways);
	
	for (i = 0; i < ninfectables; i++) si += infways[infectables[i]];

	if ((si == 0) && (nrecoverables == 0)) {
		fmpz_poly_mul(g.a, donum, wnum);
		fmpz_poly_mul(g.b, doden, wden);
		fmpz_poly_mul(g.a, g.a, g.oden);
		fmpz_poly_mul(g.onum, g.onum, g.b);
		fmpz_poly_add(g.onum, g.onum, g.a);
		fmpz_poly_mul(g.oden, g.oden, g.b);

		simplify(&g.onum, &g.oden);

		goto CLEAR_EXIT;
	}

	fmpz_poly_zero(g.b);
	fmpz_poly_set_coeff_ui(g.b, 0, (unsigned long) nrecoverables);
	fmpz_poly_set_coeff_ui(g.b, 1, (unsigned long) si);

	fmpz_poly_mul(g.a, g.b, donum);
	fmpz_poly_add(donum, g.a, doden);
	fmpz_poly_mul(doden, doden, g.b);
	fmpz_poly_mul(wden, wden, g.b);
	
	simplify(&donum, &doden);

	for (i = 0; i < ninfectables; i++) {
		you = infectables[i];
		fmpz_poly_zero(g.a);
		fmpz_poly_set_coeff_ui(g.a, 1, (unsigned long) infways[you]);
		fmpz_poly_mul(nnum, g.a, wnum);

		if (n[you].state == SENTI) {
			fmpz_poly_mul(g.a, donum, nnum);
			fmpz_poly_mul(g.b, doden, wden);
			fmpz_poly_mul(g.a, g.a, g.oden);
			fmpz_poly_mul(g.onum, g.onum, g.b);
			fmpz_poly_add(g.onum, g.onum, g.a);
			fmpz_poly_mul(g.oden, g.oden, g.b);

			simplify(&g.onum, &g.oden);

		} else {
			n[you].state = I;
			stepdown(nnum, wden, donum, doden);
			n[you].state = S;
		}
	}
	
	for (i = 0; i < nrecoverables; i++) {
		n[you = recoverables[i]].state = R;
		stepdown(wnum, wden, donum, doden);
		n[you].state = I;
	}

	CLEAR_EXIT:

	fmpz_poly_clear(donum);
	fmpz_poly_clear(doden);
	fmpz_poly_clear(wden);
	fmpz_poly_clear(nnum);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input

int main (int argc, char *argv[]) {
	int i, j;
	char s1[OSIZE], s2[OSIZE], s3[OSIZE];
	fmpz_poly_t wnum, wden, donum, doden, sonum, soden;

	if (argc < 3) {
		fprintf(stderr, "usage: ./senti [# links] [links] <seeds>\n");
		return 1;
	}

	s1[0] = s2[0] = s3[0] = '\0';
	fmpz_poly_init(g.a);
	fmpz_poly_init(g.b);
	fmpz_poly_init(wnum);
	fmpz_poly_init(wden);
	fmpz_poly_init(donum);
	fmpz_poly_init(doden);
	fmpz_poly_init(sonum);
	fmpz_poly_init(soden);
	fmpz_poly_init(g.onum);
	fmpz_poly_init(g.oden);

	// reading and setting up the network
	g.nl = atoi(argv[1]);
	g.n = 0;
	for (i = 0; i < N; i++) {
		n[i].deg = 0;
		n[i].state = S;
	}

	for (i = 0; i < g.nl; i++)
		add_edge(atoi(argv[2 + 2 * i]), atoi(argv[3 + 2 * i]));

	// setting the seeds
	for (i = 2 + 2 * g.nl, j = 0; i < argc; i++, j++) {
		n[atoi(argv[i])].state = SENTI;   
		if (j > 0) strcat(s1, " ");
		strcat(s1, argv[i]);
	}

	fmpz_poly_zero(sonum);
	fmpz_poly_one(soden);

	for (i = 0; i < g.n; i++) {
		fmpz_poly_one(wnum);
		fmpz_poly_one(wden);
		fmpz_poly_zero(donum);
		fmpz_poly_one(doden);
		fmpz_poly_zero(g.onum);
		fmpz_poly_one(g.oden);
		
		for (j = 0; j < g.n; j++) if (n[j].state != SENTI) n[j].state = S;
		if (n[i].state != SENTI) n[i].state = I;
		stepdown(wnum, wden, donum, doden);

		fmpz_poly_mul(g.a, g.onum, soden);
		fmpz_poly_mul(sonum, sonum, g.oden);
		fmpz_poly_add(sonum, sonum, g.a);
		fmpz_poly_mul(soden, soden, g.oden);

		simplify(&sonum, &soden);
	}

	fmpz_poly_scalar_mul_ui(g.b, soden, (slong) g.n);
	simplify(&sonum, &g.b);

	strc(fmpz_poly_get_str_pretty(sonum, "x"), s2);
	strc(fmpz_poly_get_str_pretty(g.b, "x"), s3);
	printf("%s, (%s)/(%s)\n", s1, s2, s3);

	fmpz_poly_clear(g.a);
	fmpz_poly_clear(g.b);
	fmpz_poly_clear(wnum);
	fmpz_poly_clear(wden);
	fmpz_poly_clear(sonum);
	fmpz_poly_clear(soden);
	fmpz_poly_clear(donum);
	fmpz_poly_clear(doden);
	fmpz_poly_clear(g.onum);
	fmpz_poly_clear(g.oden);

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
