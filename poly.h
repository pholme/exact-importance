#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fmpz_poly.h>

#define SUSCEPTIBLE 0 // flags for the node states
#define SENTINEL 1
#define INFECTIOUS 2
#define RECOVERED 3
#define VACCINATED 4
#define MAXN 8 // maximal number of nodes
#define OSIZE 10000 // for printing output

#define TRUE 1
#define FALSE 0

#define IS_SUSCEPTIBLE(x) (n[(x)].state < INFECTIOUS) // sentinels are also susceptible
#define IS_INFECTIOUS(x) (n[(x)].state == INFECTIOUS)
#define IS_RECOVERED(x) (n[(x)].state == RECOVERED)
#define IS_SENTINEL(x) (n[(x)].state == SENTINEL)
#define NOT_VACCINATED(x) (n[(x)].state != VACCINATED)
#define NOT_SENTINEL(x) (n[(x)].state != SENTINEL)

typedef struct {
	int deg, nb[MAXN], state, id;
} NODE;

typedef struct {
	int n, nn, nl, infectables[MAXN], recoverables[MAXN], check[MAXN];
	fmpz_poly_t onum, oden, a, b, s;
} GLOBAL;

void simplify (fmpz_poly_t *, fmpz_poly_t *); // in poly.c
int obsize (); // in poly.c
void add_edge (int, int); // in poly.c
void strc (char *, char *); // in poly.c
void stepdown (fmpz_poly_t, fmpz_poly_t); // in stepdown.c
void get_inf_rec (int *, int *, int *, int *, int *, int *); // in get_inf_rec.c
