#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi





//Define Initial Condition; must be funtion of variable x
double InitialCondition(double x){
    return sin(4*PI*x);
}


int main(){
	double alpha = 0.01;		//the alpha in u_t = alpha u_xx
	double x0 = 0.;			//left spatial boundary
	double xf = 1.;			//right spatial boundary
	double dt = 0.001; 	//temporal step
	int     N = 500;		//time steps
	int 	J = 500;		//space steps
	heatequation_cn(InitialCondition, alpha, x0, xf, dt, N, J);
}






