#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fmpz_poly.h>

#define N 8 // max number of nodes 
#define S 0
#define I 1
#define R 2
#define V 3
#define SENTI 4
#define OSIZE 10000

typedef struct {
	int deg, nb[N], state;
} NODE;

typedef struct {
	int n, nl;
	char s[OSIZE];
	fmpz_poly_t onum, oden, a, b;
} GLOBAL;

void simplify (fmpz_poly_t *, fmpz_poly_t *);
int obsize ();
void add_edge (int, int);
void strc (char *, char *);
