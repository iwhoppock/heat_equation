
#include "linear_algebra.h"

#include <stdlib.h>
#include <stdio.h>

// ----------------------------------------------------------------------
// matrix_create
//
// create a new matrix of size "m" x "n"

struct matrix *
matrix_create(int m, int n)
{
  struct matrix *M = malloc(sizeof(*M));
  M->m = m;
  M->n = n;
  M->vals = calloc(m * n, sizeof(*M->vals));

  return M;
}

// ----------------------------------------------------------------------
// matrix_destroy
//
// destroys (deallocates) the matrix "M"

void
matrix_destroy(struct matrix *M)
{
  free(M->vals);
  free(M);
}

// ----------------------------------------------------------------------
// matrix_print
//
// prints matrix "M" to stdout

void
matrix_print(struct matrix *M)
{
  printf("[(#=%dx%d)", M->m, M->n);
  for (int i = 0; i < M->m; i++) {
    printf(" {");
    for (int j = 0; j < M->n; j++) {
      printf(" %g", MAT(M, i, j));
    }
    printf(" }");
  }
  printf("]");
}
