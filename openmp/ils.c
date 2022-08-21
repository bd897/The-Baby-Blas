/**********************************************************************
 *
 * ITERATIVE LINEAR SOLVER
 *
 * Andrew J. Pounds, Ph.D.
 * Spring 2018
 *
 * Unless otherwise noted, all code and methods belong to the author.
 * Equations for the Jacoby iterative solver we adapted from Golub
 * and van Loan, "Matrix Computations", Johns Hopkins University press,
 * 1996. 
 *
 **********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
    void ils_( int *threads, int *len,  double *a, double *b, double *x );
#ifdef __cplusplus
}
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>
#include <omp.h>

/*  S E R I A L   C O D E  */

// Need the following prototype in case of a zero along the diagonal
void  dls_( int *threads, int *len, double *a, double *b, double *x );


void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    	int i, j, k, N, numThreads, iteration;
    	int zeros, converged;
    	double maxb, sum, sum1, sum2;
    	double ZERO = 0.0;
    	int ITERATION_MAX = 2000;
	double tol = 5.0e-15;
	numThreads = *threads;

	omp_set_num_threads(numThreads);

  	double *tmp;

    	N = *len;

	// First we are checking if there are any zeros along the diagonal
	zeros = 0;
	i=0;
	while( i<N && !zeros ){
		zeros = *(a+(i*N)+i) == ZERO;
		i++;
	}

	// If there are no zeros along the diagonal, we will use :  Iterative Jacobi Method
	// Note: We are not modifying A or B, in case we need to revert to the DLS
	if( !zeros ){
	
		// Create a temp vector to hold our values
		tmp = malloc( N * sizeof(double) );

		// Fill our temp vector with the intial values, and fill x with b
		for(i=0; i<N; i++){
			*(x+i) = 0.0;
			*(tmp+i) = *(b+i);
		}

		// If more than N/3 iterations are done, DLS is more efficient
		ITERATION_MAX = fmax(ITERATION_MAX, N/3);
			
		// Do while convergence or the iteration max has not been met
		iteration = 0;
		converged = 0;
		while( !converged && iteration < ITERATION_MAX ){
		    
			// Check for convergence, by seeing if 2-Norm is in tolerance
			maxb = *(tmp+0);
			sum = 0.0;
	 	    	for(i=0; i<N; i++){   		
				maxb = fmax(maxb, fabs(*(tmp+i)));
				sum += (*(x+i)-*(tmp+i))*(*(x+i)-*(tmp+i));
		    	}
			sum = sqrt(sum);
			converged = sum/maxb < tol;
					

		    #pragma omp parallel default(shared) private(j)
		    {
			// Copy last result to tmp
			#pragma omp for 
			for(i=0; i<N; i++) *(tmp+i) = *(x+i);
		
			#pragma omp barrier

			// Start reduction process
			#pragma omp for
			for(i=0; i<N; i++){
				sum1 = 0.0;
			
				for(j=0; j<i-1; j++) sum1 += *(a+(i*N)+j) * *(tmp+j);
				sum2 = 0.0;
				
				for(j=i+1; j<N; j++) sum2 += *(a+(i*N)+j) * *(tmp+j);
				*(x+i) = ( *(b+i) - sum1 - sum2 )/ *(a+(i*N)+i);
			}
		    }

			iteration++;
		}
		// Free our temp vector
		free(tmp);
		// If the iteration max is met, revert to the DLS
	        if ( iteration == ITERATION_MAX) {
            		printf(" *** ITERATIVE SOLVER FAILED TO REACH CONVERGENCE AFTER  ***\n");
            		printf(" *** %d ITERATIONS, SWITCHING TO DIRECT SOLVER ***\n",iteration);
            		dls_( threads, len, a, b, x );
        	}
	}
	// If there is a zero along the diagonal, revert to the DLS
	else {

        	printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        	printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        	dls_( threads, len, a, b, x );
    	}
}
