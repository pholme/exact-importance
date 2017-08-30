#include "poly.h"

extern GLOBAL g;
extern NODE n[MAXN];

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// calculated the outbreak size if no more nodes get infected

int obsize () {
	int i, s = 0;

	for (i = 0; i < g.n; i++)
		if ((IS_INFECTIOUS(i)) || (IS_RECOVERED(i))) s++;

	return s;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// building the network

void add_edge (int me, int you) {

	n[me].nb[n[me].deg++] = you;
	n[you].nb[n[you].deg++] = me;

	if (me + 1 > g.n) g.n = me + 1;
	if (you + 1 > g.n) g.n = you + 1;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// changing the FLINT style polynomial pretty print to a SymPy friendly format

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
// eliminates common factors from polys p and q, thus simplifying p/q

void simplify (fmpz_poly_t *p, fmpz_poly_t *q) {

	fmpz_poly_gcd(g.s, *p, *q);
	fmpz_poly_div(*p, *p, g.s);
	fmpz_poly_div(*q, *q, g.s);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
