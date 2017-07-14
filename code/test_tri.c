#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi




int main(){
    struct matrix *A = matrix_create(3,3);
    MAT(A,0,0) = 1.;
    MAT(A,1,1) = 2.;
    MAT(A,2,2) = 3.;
    MAT(A,0,1) = 4.;
    MAT(A,0,2) = 0.;
    MAT(A,1,2) = 6.;
    MAT(A,1,0) = 7.;
    MAT(A,2,0) = 0.;
    MAT(A,2,1) = 9.;
    //matrix_print(A);
    //
    //  1   4   0                   1       0
    //  7   2   6   (inverse)  *    2   =   0.25
    //  0   9   3                   3       0.25
    //
    struct vector *RHS =vector_create_and_set(3, (double[3]) {1.,2.,3.});
    struct vector *a =vector_create_and_set(3, (double[3]) {0.,7.,9.});
    struct vector *b =vector_create_and_set(3, (double[3]) {1.,2.,3.});
    struct vector *c =vector_create_and_set(3, (double[3]) {4.,6.,0.});
    tri(RHS,a,b,c);
    /*
    for(int i=0; i<3; i++){
        printf("%g\n",RHS->vals[i]);
    }
    */

    assert(RHS->vals[0]==0.);
    assert(RHS->vals[1]==0.25);
    assert(RHS->vals[2]==0.25);

}

















