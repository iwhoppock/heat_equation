#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi


//funtion for two norm
double norm(struct vector *x){  
	double m_sum = 0.0;
	for (int i = 0; i < x->n; ++i){
    	m_sum += x->vals[i] * x->vals[i];
	}
	return sqrt(m_sum);
}