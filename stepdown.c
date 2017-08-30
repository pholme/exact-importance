// This file includes the main function stepping down the tree of
// configurations. It is called by infmax and vacc. senti needs a separate
// and more elaborate function.

#include "poly.h"

extern GLOBAL g;
extern NODE n[MAXN];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void stepdown (fmpz_poly_t wnum, fmpz_poly_t wden) {
	int i, you, multiplicity[MAXN], to_step_to[MAXN], ninfectables;
	int nto_step_to, nr, nsi;
	fmpz_poly_t a, b, den;

	// get the configurations to go down to in the next step
	// isomorphic configurations that are possible to reach in many ways
	// are only needs to be evaluated once (then their probabilities
	// multiplied by multiplicity)

	get_inf_rec(&ninfectables, to_step_to, &nto_step_to, multiplicity, &nr,
			&nsi);

	if (nsi == 0) { // if there are no SI edges we can quit, the outbreak
		// size will be the nuumber of I + R
		fmpz_poly_mul(g.onum, g.onum, wden);
		fmpz_poly_mul(g.a, g.oden, wnum);
		fmpz_poly_scalar_mul_ui(g.a, g.a, (unsigned long) obsize());
		fmpz_poly_add(g.onum, g.onum, g.a);
		simplify(&g.onum, &g.oden);
		fmpz_poly_gcd(g.s, g.onum, wden);
		fmpz_poly_div(g.a, wden, g.s);
		fmpz_poly_div(g.onum, g.onum, g.s);
		fmpz_poly_mul(g.oden, g.oden, g.a);

		return;
	}

	// calculate the new denominator for the probability of the
	// configurations
	fmpz_poly_init(a);
	fmpz_poly_init(b);
	fmpz_poly_init(den);
	fmpz_poly_set_ui(den, (unsigned long) nr);
	fmpz_poly_set_coeff_ui(den, 1, (unsigned long) nsi);

	fmpz_poly_mul(b, wden, den);

	// go over the infection events
	for (i = 0; i < ninfectables; i++) {
		you = to_step_to[i];

		n[you].state = INFECTIOUS;

		fmpz_poly_zero(a);
		fmpz_poly_set_coeff_ui(a, 1, (unsigned long) multiplicity[you]);
		fmpz_poly_mul(a, a, wnum);

		stepdown(a, b);

		n[you].state = SUSCEPTIBLE;
	}
	
	// go over the recovery events
	for ( ; i < nto_step_to; i++) {
		you = to_step_to[i];

		n[you].state = RECOVERED;

		fmpz_poly_scalar_mul_ui(a, wnum, (unsigned long) multiplicity[you]);

		stepdown(a, b);

		n[you].state = INFECTIOUS;
	}

	fmpz_poly_clear(den);
	fmpz_poly_clear(a);
	fmpz_poly_clear(b);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
