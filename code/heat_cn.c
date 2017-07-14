#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linear_algebra.h"

#define PI 4*atan(1) //a precise calculation of Pi



void heatequation_cn(double InitialCondition(double x), double alpha, double x0, double xf, double dt, int N, int J){
	int i;
	double dx = (xf - x0) / (J - 1.); //spatial step
	struct vector *x = vector_create(J);
	for (i = 0; i < x->n; i++){
		x->vals[i] = i * dx;
	}
	struct vector *t = vector_create(N+1);
	for(i = 0; i < N+1; i++){ //N+1 because we are adding one more below
		t->vals[i] = i * dt;
	}

	////This prints spacesteps to be used for x-axis plotting////
		FILE *fout1 = NULL;									 ////
		fout1 = fopen("x-axis_spacesteps.txt", "w");		 ////
		for(i=0; i<x->n;i++){								 ////
			fprintf(fout1, "%.17lg\n", x->vals[i]);			 ////
		}													 ////
		fclose(fout1);										 ////
	////This prints timesteps to be used for x-axis plotting ////
		FILE *fout4 = NULL;									 ////
		fout4 = fopen("x-axis_timesteps.txt", "w");			 ////
		for(i=0; i<t->n;i++){								 ////
			fprintf(fout4, "%.17lg\n", t->vals[i]);			 ////
		}													 ////
		fclose(fout4);										 ////
	/////////////////////////////////////////////////////////////

	

	//initial condition @t=0
	struct vector *ic = vector_create(J-1);
	for (i = 0; i < J; i++){
		ic->vals[i] = InitialCondition(x->vals[i]);	
	}
	
	/*
	//CHECK ALL OF THE ABOVE: (prints all the above)
	for (i = 0; i < J; i++){
		printf("%g\t%g\t%g\n", x->vals[i], ic->vals[i],t->vals[i]);
	}
	*/	
	

	//define ratio r to make calculations a bit easier
	double r = alpha * dt / (2 * dx * dx);
	//CHECK:
	//printf("%g\n",r);




	//Since we will be using a tridiagonal matrix for all applications, 
	//we need only three vectors: create vectors a, b, and c. 
	//These are the sub, main and superdiagonals, respectively, of 
	//our matrix. See tridiagonal matrix function for explaination 
	//of indexing
	
	struct vector *a = vector_create(J);
	a->vals[0]=0.;
	for(i=1;i < J; i++){
		a->vals[i] = -r; 
	}
	
	struct vector *b = vector_create(J); 
	for(i=0;i<J;i++){
		b->vals[i] = 1 + 2 * r; 
	} 
	
	struct vector *c = vector_create(J);	
	for(i=0;i<J-1;i++){
		c->vals[i] = -r;
	} 
	c->vals[J-1]=0.;
	
	/*
	//CHECK
	for(i=0; i<J; i++){
		printf("%g\t%g\t%g\n",a->vals[i],b->vals[i],c->vals[i]);
	}
	*/
	


	//CRANK-NICOLSON
	//appropriately named right hand side of Ax=b problem
	struct vector *RHS = vector_create(J);
	RHS->vals[0] = 0.; 
	for(i = 1; i < J-1; i++){
		RHS->vals[i] = ic->vals[i] + r * (ic->vals[i+1] - 2 * ic->vals[i] + ic->vals[i-1]);
	}
	RHS->vals[J-1] = 0.;

	//CHECK
	/*
	for(i=0; i<J; i++){
	printf("%g\n", RHS->vals[i]);
	}
	*/
	


	//OUTPUT FILE
	FILE *fout2 = NULL;
	fout2 = fopen("results_diffusion.txt", "w");

	FILE *fout3 = NULL;
	fout3 = fopen("results_norms.txt", "w");

	//print initial state:
	for(i=0; i<ic->n;i++){
		fprintf(fout2, "%.17le\t", ic->vals[i]);
	}
	fprintf(fout2, "\n");

	fprintf(fout3, "%g\n",norm(ic));

	//TIME LOOP
	for(int j=0; j < N; j++){
		//tridiagonal algorithm. Much faster than LU 
		tri(RHS,a,b,c);
		RHS->vals[J-1] = 0.;
		RHS->vals[0]=0.;
		/*
		//Print just one iteration of RHS (must comment out time loop)
		for(i=0; i<RHS->n-1; i++){
			printf("%g\n",RHS->vals[i]);
		}
		*/
		//PRINTS SOLUTION IN ROWS. THE USER MUST PLOT BY ROWS OR SIMPLY TAKE TRANSPOSE
		for(i=0; i < RHS->n-1; i++){
			fprintf(fout2, "%.17le\t", RHS->vals[i]);
		}
		fprintf(fout2, "\n"); //next row

		//prints out norm of each solution
		fprintf(fout3, "%e\n", norm(RHS));
		

	}
	fclose(fout2);
	fclose(fout3);


	vector_destroy(x);
    x = NULL;
    vector_destroy(ic);
    ic = NULL;
    vector_destroy(a);
    a = NULL;
    vector_destroy(b);
    b = NULL;
    vector_destroy(c);
    c = NULL;
    vector_destroy(RHS);
    RHS = NULL;
}