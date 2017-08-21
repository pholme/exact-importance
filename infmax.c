// C code corresponding to infmax.py (using the FLINT library)
// see infmax.py for more detailed comments

#include "poly.h"

NODE n[N];
GLOBAL g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// making a list of nodes that can be infected or recovered
// (very naive, not the bottleneck anyway)

void get_infect_reco (int *infectables, int *ninfectables,
		int *recoverables, int *nrecoverables, int *infways) {
	int i, me, you;

	for (i = 0; i < g.n; i++) infways[i] = 0;

	for (me = 0; me < g.n; me++) {
		if (n[me].state == S) {
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

void stepdown (fmpz_poly_t wnum, fmpz_poly_t wden) {
	int i, si = 0, you, infways[N];
	int infectables[N], ninfectables = 0, recoverables[N], nrecoverables = 0;
	fmpz_poly_t a, b, den;

	get_infect_reco(infectables, &ninfectables, recoverables, &nrecoverables, infways);
	
	for (i = 0; i < ninfectables; i++) si += infways[infectables[i]];

	if (si == 0) {
		fmpz_poly_mul(g.onum, g.onum, wden);
		fmpz_poly_mul(g.a, g.oden, wnum);
		fmpz_poly_scalar_mul_ui(g.a, g.a, (slong) obsize());
		fmpz_poly_add(g.onum, g.onum, g.a);
		fmpz_poly_mul(g.oden, g.oden, wden);

		simplify(&g.onum, &g.oden);

		return;
	}

	fmpz_poly_init(a);
	fmpz_poly_init(b);
	fmpz_poly_init(den);
	fmpz_poly_set_ui(den, (unsigned long) nrecoverables);
	fmpz_poly_set_coeff_ui(den, 1, (unsigned long) si);

	for (i = 0; i < ninfectables; i++) {
		you = infectables[i];

		n[you].state = I;

		fmpz_poly_zero(a);
		fmpz_poly_set_coeff_ui(a, 1, (unsigned long) infways[you]);
		fmpz_poly_mul(a, a, wnum);
		fmpz_poly_mul(b, wden, den);

		stepdown(a, b);

		n[you].state = S;
	}
	
	for (i = 0; i < nrecoverables; i++) {
		you = recoverables[i];

		n[you].state = R;

		fmpz_poly_mul(a, wden, den);

		stepdown(wnum, a);

		n[you].state = I;
	}

	fmpz_poly_clear(den);
	fmpz_poly_clear(a);
	fmpz_poly_clear(b);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input

int main (int argc, char *argv[]) {
	int i, j;
	char s1[OSIZE], s2[OSIZE], s3[OSIZE];
	fmpz_poly_t wnum, wden;

	if (argc < 3) {
		fprintf(stderr, "usage: ./infmax [# links] [links] <seeds>\n");
		return 1;
	}

	s1[0] = s2[0] = s3[0] = '\0';
	fmpz_poly_init(wnum);
	fmpz_poly_init(wden);
	fmpz_poly_init(g.onum);
	fmpz_poly_init(g.oden);
	fmpz_poly_init(g.a);

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
		n[atoi(argv[i])].state = I;   
		if (j > 0) strcat(s1, " ");
		strcat(s1, argv[i]);
	}

	fmpz_poly_one(wnum);
	fmpz_poly_one(wden);
	fmpz_poly_zero(g.onum);
	fmpz_poly_one(g.oden);

	stepdown(wnum, wden);
	simplify(&g.onum, &g.oden);
	
	strc(fmpz_poly_get_str_pretty(g.onum, "x"), s2);
	strc(fmpz_poly_get_str_pretty(g.oden, "x"), s3);
	printf("%s, (%s)/(%s)\n", s1, s2, s3);

	fmpz_poly_clear(wnum);
	fmpz_poly_clear(wden);
	fmpz_poly_clear(g.onum);
	fmpz_poly_clear(g.oden);
	fmpz_poly_clear(g.a);

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
