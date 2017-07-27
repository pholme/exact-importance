#include "poly.h"

extern NODE n[N];
extern GLOBAL g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int obsize () {
	int i, si = 0;

	for (i = 0; i < g.n; i++) if ((n[i].state == I) || (n[i].state == R)) si++;

	return si;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void add_edge (int me, int you) {

	n[me].nb[n[me].deg++] = you;
	n[you].nb[n[you].deg++] = me;

	if (me + 1 > g.n) g.n = me + 1;
	if (you + 1 > g.n) g.n = you + 1;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void strc (char *sin, char *sout) {
	int l;
	char *ch = sin;

	while ((*ch != '\0') && ((*ch != '\n'))) {
		if (*ch == '^') strcat(sout, "**");
		else {
			l = strlen(sout);
			sout[l] = *ch;
			sout[l + 1] = '\0';
		}
		ch++;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void simplify (fmpz_poly_t *p, fmpz_poly_t *q) {
	fmpz_poly_t p0, q0, s;

	fmpz_poly_init(p0);
	fmpz_poly_init(q0);
	fmpz_poly_init(s);

	fmpz_poly_set(p0, *p);
	fmpz_poly_set(q0, *q);

	fmpz_poly_gcd_modular(s, p0, q0);

	fmpz_poly_div(*p, p0, s);
	fmpz_poly_div(*q, q0, s);

	fmpz_poly_clear(p0);
	fmpz_poly_clear(q0);
	fmpz_poly_clear(s);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -