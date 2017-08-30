// C code corresponding to senti.py (using the FLINT library)

#include "poly.h"
#include <fmpq.h>

GLOBAL g;
NODE n[MAXN];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void stepdown_senti (fmpz_poly_t wnum, fmpz_poly_t wden, fmpz_poly_t donum,
		fmpz_poly_t doden) {
	int i, you, multiplicity[MAXN], to_step_to[MAXN], ninfectables;
	int nto_step_to, nr, nsi;
	fmpz_poly_t donum0, doden0, wden0, wnum0;
	fmpq_t q;

	get_inf_rec(&ninfectables, to_step_to, &nto_step_to, multiplicity, &nr,
			&nsi);

	if (nsi == 0) { // if no SI links the outbreak will die
		if (nr > 0) {
			fmpq_init(q);
			fmpq_harmonic_ui(q, nr); // time to extinction from here is a partial harmonic sum
			fmpz_poly_scalar_mul_fmpz(g.a, doden, fmpq_numref(q));
			fmpz_poly_scalar_mul_fmpz(g.b, donum, fmpq_denref(q));
			fmpz_poly_add(g.a, g.a, g.b);
			fmpz_poly_scalar_mul_fmpz(g.b, doden, fmpq_denref(q));
			fmpq_clear(q);
			simplify(&g.a, &g.b); // faster to simplify first, multiply later
			fmpz_poly_mul(g.a, g.a, wnum);
			fmpz_poly_mul(g.b, g.b, wden);
		} else {
			fmpz_poly_mul(g.a, donum, wnum);
			fmpz_poly_mul(g.b, doden, wden);
		}

		fmpz_poly_mul(g.a, g.a, g.oden);
		fmpz_poly_mul(g.onum, g.onum, g.b);
		fmpz_poly_add(g.onum, g.onum, g.a);

		simplify(&g.onum, &g.oden);
		simplify(&g.onum, &g.b);
		fmpz_poly_mul(g.oden, g.oden, g.b);

		return;
	}
	
	// calculating denominators
	fmpz_poly_init2(donum0, 1 + fmpz_poly_length(donum));
	fmpz_poly_init2(doden0, 1 + fmpz_poly_length(doden));
	fmpz_poly_init2(wnum0, 1 + fmpz_poly_length(wnum));
	fmpz_poly_init2(wden0, 1 + fmpz_poly_length(wden));

	fmpz_poly_set_ui(g.a, (unsigned long) nr);
	fmpz_poly_set_coeff_ui(g.a, 1, (unsigned long) nsi);

	fmpz_poly_mul(wden0, wden, g.a);

	fmpz_poly_mul(donum0, g.a, donum);
	fmpz_poly_add(donum0, donum0, doden);

	fmpz_poly_gcd(g.s, donum0, doden);
	fmpz_poly_div(donum0, donum0, g.s);
	fmpz_poly_div(g.b, doden, g.s);

	fmpz_poly_gcd(g.s, donum0, g.a);
	fmpz_poly_div(donum0, donum0, g.s);
	fmpz_poly_div(g.a, g.a, g.s);

	fmpz_poly_mul(doden0, g.b, g.a);

	for (i = 0; i < ninfectables; i++) { // going over the infection events
		you = to_step_to[i];
		fmpz_poly_zero(g.a);
		fmpz_poly_set_coeff_ui(g.a, 1, (unsigned long) multiplicity[you]);
		fmpz_poly_mul(wnum0, g.a, wnum);

		if (IS_SENTINEL(you)) {
			fmpz_poly_mul(g.a, donum0, wnum0);
			fmpz_poly_mul(g.b, doden0, wden0);
			fmpz_poly_mul(g.a, g.a, g.oden);
			fmpz_poly_mul(g.onum, g.onum, g.b);
			fmpz_poly_add(g.onum, g.onum, g.a);

			simplify(&g.onum, &g.oden);
			simplify(&g.onum, &g.b);

			fmpz_poly_mul(g.oden, g.oden, g.b);
		} else {
			n[you].state = INFECTIOUS;
			stepdown_senti(wnum0, wden0, donum0, doden0);
			n[you].state = SUSCEPTIBLE;
		}
	}
	
	for ( ; i < nto_step_to; i++) {  // going over the recovery events
		you = to_step_to[i];

		n[you].state = RECOVERED;
		fmpz_poly_scalar_mul_ui(wnum0, wnum, (unsigned long) multiplicity[you]);
		stepdown_senti(wnum0, wden0, donum0, doden0);
		n[you].state = INFECTIOUS;
	}

	fmpz_poly_clear(donum0);
	fmpz_poly_clear(doden0);
	fmpz_poly_clear(wden0);
	fmpz_poly_clear(wnum0);
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
	fmpz_poly_init(g.s);
	// (wnum / wden) is the probability of a branch of the search tree 
	fmpz_poly_init(wnum);
	fmpz_poly_init(wden);
	// (donum / doden) is the duration to traverse a branch of the search tree 
	fmpz_poly_init(donum);
	fmpz_poly_init(doden);
	// (g.onum / g.oden) is the weighted sum of durations over branches
	fmpz_poly_init(g.onum);
	fmpz_poly_init(g.oden);
	// (sonum / soden) is the sum over seeds of (g.onum / g.oden)
	fmpz_poly_init(sonum);
	fmpz_poly_init(soden);

	// reading and setting up the network
	g.nl = atoi(argv[1]);
	g.n = 0;
	for (i = 0; i < MAXN; i++) {
		n[i].deg = 0;
		n[i].id = i;
		n[i].state = SUSCEPTIBLE;
	}

	for (i = 0; i < g.nl; i++)
		add_edge(atoi(argv[2 + 2 * i]), atoi(argv[3 + 2 * i]));

	// setting the seeds
	for (i = 2 + 2 * g.nl, j = 0; i < argc; i++, j++) {
		n[atoi(argv[i])].state = SENTINEL;   
		if (j > 0) strcat(s1, " ");
		strcat(s1, argv[i]);
	}

	fmpz_poly_zero(sonum);
	fmpz_poly_one(soden);

	for (i = 0; i < g.n; i++) { // loop over all nodes as seeds
		fmpz_poly_one(wnum);
		fmpz_poly_one(wden);
		fmpz_poly_zero(donum);
		fmpz_poly_one(doden);
		fmpz_poly_zero(g.onum);
		fmpz_poly_one(g.oden);
		
		for (j = 0; j < g.n; j++) if (NOT_SENTINEL(j)) n[j].state = SUSCEPTIBLE;
		if (NOT_SENTINEL(i)) n[i].state = INFECTIOUS;
		stepdown_senti(wnum, wden, donum, doden);

		fmpz_poly_mul(g.a, g.onum, soden);
		fmpz_poly_mul(sonum, sonum, g.oden);
		fmpz_poly_add(sonum, sonum, g.a);
		simplify(&sonum, &soden);
		simplify(&sonum, &g.oden);
		fmpz_poly_mul(soden, soden, g.oden);
	}

	fmpz_poly_scalar_mul_ui(g.b, soden, (slong) g.n);
	simplify(&sonum, &g.b);

	strc(fmpz_poly_get_str_pretty(sonum, "x"), s2);
	strc(fmpz_poly_get_str_pretty(g.b, "x"), s3);
	printf("%s, (%s)/(%s)\n", s1, s2, s3);

	fmpz_poly_clear(g.a);
	fmpz_poly_clear(g.b);
	fmpz_poly_clear(g.s);
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
