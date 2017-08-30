// calculating the expected outbreak size after some nodes are vaccinated

#include "poly.h"

NODE n[MAXN];
GLOBAL g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input

int main (int argc, char *argv[]) {
	int i, j;
	char s1[OSIZE], s2[OSIZE], s3[OSIZE];
	fmpz_poly_t wnum, wden, sonum, soden;

	if (argc < 3) {
		fprintf(stderr, "usage: ./vacc [# links] [links] <seeds>\n");
		return 1;
	}

	s1[0] = s2[0] = s3[0] = '\0';
	fmpz_poly_init(wnum);
	fmpz_poly_init(wden);
	fmpz_poly_init(g.onum);
	fmpz_poly_init(g.oden);
	fmpz_poly_init(sonum);
	fmpz_poly_init(soden);
	fmpz_poly_init(g.a);
	fmpz_poly_init(g.s);

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
		n[atoi(argv[i])].state = VACCINATED;   
		if (j > 0) strcat(s1, " ");
		strcat(s1, argv[i]);
	}

	fmpz_poly_zero(sonum);
	fmpz_poly_one(soden);

	for (i = 0; i < g.n; i++) if (NOT_VACCINATED(i)) {
		fmpz_poly_one(wnum);
		fmpz_poly_one(wden);
		fmpz_poly_zero(g.onum);
		fmpz_poly_one(g.oden);

		for (j = 0; j < g.n; j++)
			if (NOT_VACCINATED(j)) n[j].state = SUSCEPTIBLE;
		n[i].state = INFECTIOUS;
		stepdown(wnum, wden);

		fmpz_poly_mul(g.a, g.onum, soden);
		fmpz_poly_mul(sonum, sonum, g.oden);
		fmpz_poly_add(sonum, sonum, g.a);
		simplify(&sonum, &soden);
		simplify(&sonum, &g.oden);
		fmpz_poly_mul(soden, soden, g.oden);
	}

	fmpz_poly_scalar_mul_ui(g.a, soden, (slong) g.n);
	simplify(&sonum, &g.a);

	strc(fmpz_poly_get_str_pretty(sonum, "x"), s2);
	strc(fmpz_poly_get_str_pretty(g.a, "x"), s3);
	printf("%s, (%s)/(%s)\n", s1, s2, s3);

	fmpz_poly_clear(wnum);
	fmpz_poly_clear(wden);
	fmpz_poly_clear(g.onum);
	fmpz_poly_clear(g.oden);
	fmpz_poly_clear(sonum);
	fmpz_poly_clear(soden);
	fmpz_poly_clear(g.a);
	fmpz_poly_clear(g.s);

	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
