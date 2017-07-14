#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi




//standard tridiagonal algorithm taken from open source and modifed to fit my macros
void tri(struct vector *x, struct vector *a, struct vector *b, struct vector *c){
 	 /*
     solves Ax = v where A is a tridiagonal matrix consisting a, b, c
     x - Input vector v, returns the solution x
     a - subdiagonal: <0, mat(1,0),..., mat(n+1, n)>
     b - the main diagonal: <mat(1,1),..., mat(n,n)>
     c - superdiagonal: <mat(0,1), ..., mat(n,n+1), 0>
     */
    c->vals[0] = c->vals[0] / b->vals[0];
    x->vals[0] = x->vals[0] / b->vals[0];
    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (int ix = 1; ix < b->n; ix++) {
        const double m = 1. / (b->vals[ix] - a->vals[ix] * c->vals[ix - 1]);
        c->vals[ix] = c->vals[ix] * m;
        x->vals[ix] = (x->vals[ix] - a->vals[ix] * x->vals[ix - 1]) * m;
    }   
    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (int ix = b->n - 2; ix >= 0; ix--){
        x->vals[ix] -= c->vals[ix] * x->vals[ix + 1];
    }
}