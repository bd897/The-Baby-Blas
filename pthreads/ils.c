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

#include <pthread.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>

/* Pthreads ILS */

// Need the following prototype in case of a zero along the diagonal
void  dls_( int *threads, int *len, double *a, double *b, double *x );
void  *ils_thread_worker();

struct args{
	int N;
	int startRow;
	int stopRow;
	double *Aptr;
	double *Bptr;
	double *Cptr;
	double *Dptr;
};

void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    	int i, j, k, N, iteration, startRow, stopRow;
    	int zeros, converged;
    	double maxb, sum, sum1, sum2;
    	double ZERO = 0.0;
	int numThreads = *threads;
    	int ITERATION_MAX = 2000;
	double tol = 5.0e-15;
	struct args *thread_args;
  	double *tmp;
	int  *numRows;
	pthread_t *thread_id;
    	N = *len;
	
	// Allocate memory for pthread variables
	thread_id = (pthread_t *) malloc (numThreads*sizeof(pthread_t));
	numRows = (int *)malloc(numThreads*sizeof(int));

	// Divide up the rows to be worked
	for(i=0; i<numThreads;i++){
		*(numRows+i) = N/numThreads;	
	}
	for(i=0; i<N%numThreads; i++){
		*(numRows+i) = *(numRows+i) + 1;
	}

	// We are checking if there are any zeros along the diagonal
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

			// Copy last result to tmp
			for(i=0; i<N; i++) *(tmp+i) = *(x+i);
			
			// Make thread args and send off threads
			stopRow = 0;
			for(i=0;i<numThreads; i++){
			    {
				startRow = stopRow;
				stopRow = startRow + *(numRows+i);
				thread_args = (struct args *)malloc(sizeof(struct args));
				thread_args -> N = N;
				thread_args -> startRow = startRow;
				thread_args -> stopRow = stopRow;
				thread_args -> Aptr = a;
				thread_args -> Bptr = b;
				thread_args -> Cptr = tmp;
				thread_args -> Dptr = x;
				pthread_create( thread_id+i, NULL, &ils_thread_worker, thread_args);
			    }
			}
			// Wait for all threads, then move on
			for(i=0;i<numThreads;i++){
				pthread_join( *(thread_id+i), NULL);
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
		printf(" Iteration Count : %d \n", iteration);
	}
	// If there is a zero along the diagonal, revert to the DLS
	else {

        	printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        	printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        	dls_( threads, len, a, b, x );
	}
	free(numRows);
	free(thread_id);    	
}

void *ils_thread_worker( struct args *thread_args){


	int i,j;
	double sum1, sum2;
	int NN, rowStart, rowStop, kk;
	double *a, *b, *tmp, *x;	

	// Get Struct Args
	NN = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop = thread_args->stopRow;
	a = thread_args->Aptr;
	b = thread_args->Bptr;
	tmp = thread_args->Cptr;
	x = thread_args->Dptr;

	// Start reduction process
	for(i=rowStart; i<rowStop; i++){
		sum1 = 0.0;
		for(j=0; j<i-1; j++) sum1 += *(a+(i*NN)+j) * *(tmp+j);
        	sum2 = 0.0;
        	for(j=i+1; j<NN; j++) sum2 += *(a+(i*NN)+j) * *(tmp+j);
        	*(x+i) = ( *(b+i) - sum1 - sum2 )/ *(a+(i*NN)+i);
	}
	// Clean up
	free(thread_args);
	// Send thread back
	pthread_exit(NULL);
}
