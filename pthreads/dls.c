#ifdef __cplusplus
extern "C" {
#endif
    void dls_( int *threads, int *len,  double *A, double *b, double*rvec);
#ifdef __cplusplus
}
#endif

#include <pthread.h>
#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>

// Prototype
void *dls_thread_worker();

// Struct for threaded reduction
struct args{
	int N;
	int startRow;
	int stopRow;
	int k;
	double *Aptr;
};


// A direct linear solver to find the vector x in Ax=b
void dls_(int *threads, int *len, double *A, double *b, double *rvec){

	int i, j, k, u, rows, rows2, N;
	double sum, pivotMax, tmp, ZERO;
	int iPivot, numThreads, startRow, stopRow;
	double *v;
	int *p, *numRows;
	int isSymmetric, posDef, zeros;
	pthread_t *thread_id;
	struct args *thread_args;	

	// Get N size
	N = *len;
	// Get number of threads
	numThreads = *threads;
	// Initialize variables
	ZERO = 0.0;
	isSymmetric = 1;		
	posDef = 1;
	zeros = 1;

	// Checks if our matrix is symmetric
	while(i<N && isSymmetric){
		for(i=0; i<N; i++){
			for(j=0; j<i; j++){
				isSymmetric = *(A+(i*N)+j) == *(A+(j*N)+i);
			}
		}
	}
	// If our matrix is not symmetric, do PLU Factorization
	if( !isSymmetric ){
	
		// Allocate memory for pivot array and fill it
		p = malloc( (N-1) * sizeof(int) );
		for (k=0;k<N-1;k++) *(p+k)=k;       

		// Find the greatest value in the column and record its position
		for (k=0;k<N-1;k++) {
			pivotMax = *(A+(k*N)+k);
		        iPivot = k; 
		        for (u=k;u<N;u++) {
		        	if ( fabs(*(A+u*N+k)) > fabs(pivotMax) ) {
		                	pivotMax = *(A+u*N+k);
					iPivot = u;
		                }        
			}

			// If a greater pivot value is found, swap rows, and record pivot
			if ( iPivot != k ) {
		        	u = iPivot; 
		                for (j=k;j<N;j++) {
		                      	tmp = *(A+k*N+j);
		                       	*(A+(k*N)+j) = *(A+u*N+j);
		         		*(A+(u*N)+j) = tmp;
				}
			} 
			*(p+k) = iPivot;

			// If the pivot is not zero, reduce and eliminate the column
			if( *(A+(k*N)+k) != ZERO ) {
			
				// Divide up work for threads			    
			    	if( (N-(k+1)) < numThreads ){
			    		for (rows=k+1;rows<N;rows++) { 
                        			*(A+(rows*N)+k) = *(A+(rows*N)+k) / *(A+(k*N)+k);
						for (rows2=k+1;rows2<N;rows2++) { 
                        	    			*(A+(rows*N)+rows2) = *(A+(rows*N)+rows2) - *(A+(rows*N)+k) * *(A+(k*N)+rows2) ;
                    	    			}
                	    		}
				}
				else{
					// Allocate memory for thread arg
					thread_id = (pthread_t *) malloc (numThreads*sizeof(pthread_t));
					numRows = (int *)malloc(numThreads*sizeof(int));
					// Divide up rows to be worked
					for(i=0; i<numThreads;i++){
						*(numRows+i) = (N-(k+1))/numThreads;	
					}
					for(i=0; i< (N-(k+1))%numThreads; i++){
						*(numRows+i) = *(numRows+i) + 1;
					}
					// Assign values to struct and send off threads
					stopRow = k+1;
					for(i=0;i<numThreads; i++){
					    {
						startRow = stopRow;
						stopRow = startRow + *(numRows+i);
						thread_args = (struct args *)malloc(sizeof(struct args));
						thread_args -> N = N;
						thread_args -> startRow = startRow;
						thread_args -> stopRow = stopRow;
						thread_args -> k = k;
						thread_args -> Aptr = A;
						pthread_create( thread_id+i, NULL, &dls_thread_worker, thread_args);
					    }
					}
					// Wait for all threads to be done before moving on
					for(i=0;i<numThreads;i++){
						pthread_join( *(thread_id+i), NULL);
					}
			        }
			}
			// If else, the matrix is singular and stop the program
			else{
		    	    printf( "Element a[%d][%d} = %f\n", k, k, *(A+k*N+k)); 
                    	    printf( " *** MATRIX A IS SINGULAR *** \n");
	                    printf( "    -- EXECUTION HALTED --\n");
        	            exit(1);
			}
	 	}
		// Swaps rows of b using p, to apply the same operations as A
		for(k=0; k<N-1; k++){
			// Swap rows of x with p(k)
			tmp = *(b+k);
            		*(b+k) = *(b+ *(p+k));
	                *(b+ *(p+k)) = tmp;
			// Solves Ly = b
            		for (j=k+1;j<N;j++){
                		*(b+j)= *(b+j) - *(b+k) * *(A+N*j+k);  
			}
		}
		// Now we are in Ux = y form and can do back sub. to solve for x
        	*(b+N-1) = *(b+N-1) / *(A+N*(N-1)+(N-1));
        	for (i=N-2;i>=0;i--){
            		tmp = 0.0;
            		for (j=i+1;j<N;j++) {
                		tmp = tmp + *(A+i*N+j) * *(b+j);
            		}
           		*(b+i) = ( *(b+i) - tmp ) / *(A+i*N+i); 
        	}
		// Fill our results vector with the values in b
		for (i=0;i<N;i++){
			*(rvec+i) = *(b+i);
		}
		free(p);
	}	
	// If symmetric, do an LDL^t factorization
	else{	
		// Checks if our matrix is singular, by checking diagonal for zeros
		zeros = 1;
		i = 0;
		while( i < N && zeros) {
			zeros = *(A+(i*N)+i) == ZERO;
		}

		// If singular, the program is stopped
		if( zeros ){
				
            		printf( " *** MATRIX A IS SINGULAR *** \n");
            		printf( "    -- EXECUTION HALTED -- \n");
            		exit(1);		
		}
		// If not singular, continue with LDL^t
		// Note: We are storing L in A and D in rvec during factorization
		else{
			// Allocate memory to V
	                v = malloc( (N) * sizeof(double) );
			
			// While doing LDL^t, we can see if positive definite
			// by checking the diagonal for any negative entries
			
			// We are storing D in rvec and L in A
			for(i=0; i<N; i++){
				sum = 0.0;
				for(j=0; j<i; j++){
					*(v+j) = *(A+(i*N)+j) * *(rvec+j);
					sum += *(A+(i*N)+j) * *(v+j);
				}

				// If a diagonal entry is less than zero, not pos def	
				*(rvec+i) = *(A+(i*N)+i);
				if(posDef){
					posDef = *(rvec+i) >= ZERO;
				}
				// Now we are solving for L
				for(j=i+1; j<N; j++){
					sum = 0.0;
					for(k=0; k<i; k++){
						sum += *(A+(j*N)+k) * *(v+k);
					}
					*(A+(j*N)+i) = ( *(A+(i*N)+j) - sum)/ *(rvec+i);
				}
			}
	
			// If not positive definite solve from LDL^t
			if( !posDef ){

				// Solve for Lz = b, store z in v vector from earlier
				for(i=0; i<N; i++){
					sum = 0.0;
					for(j=0; j<i; j++){
						sum += *(A+(i*N)+j) * *(v+j);
					}	
					*(v+i) = *(b+i) - sum;
				}				

				// Solve for Dy = z, we will store y in D
				for(i=0; i<N; i++){
					*(rvec+i) = *(v+i) / *(rvec+i);
				}

				// Solve for L^t x = y
				sum = 0.0;
				for(i=0; i<N; i++){
					for(j=i+1; j<N; j++){
						sum += *(A+(i*N)+k) * *(rvec+j);
					}
					*(rvec+i) = *(rvec+i) - sum;
				}
				// Our answer should be in rvec now
			}
			// Now we know the matrix is pos def, and we can do Cholesky
			else{
				// Solve for Uy = b
				for(i=0; i<N; i++){
					sum = 0.0;
					for(j=0; i<i; j++){
						sum += *(A+(i*N)+j) * *(rvec+j);
					}
					*(rvec+i) = (*(b+i) - sum)/ *(A+(i*N)+i);
				}		
				// Solve for U^t x = y
				for(i=0; i<N; i++){
					sum = 0.0;
					for(j=i+1; j<N; j++){
						sum += *(A+(j*N)+i) * *(rvec+j);
					}	
					*(rvec+i) = ( *(rvec+i) - sum ) / *(A+(i*N)+i);
				}					
			}
		}
		free(v);
	}
}



void *dls_thread_worker( struct args *thread_args){

	int rows, rows2;
	int NN, rowStart, rowStop, kk;
	double *mat;	

	// Get Struct Args
	NN = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop = thread_args->stopRow;
	kk = thread_args->k;
	mat = thread_args->Aptr;

	for (rows=rowStart;rows<rowStop;rows++) { 
		*(mat+(rows*NN)+kk) = *(mat+(rows*NN)+kk) / *(mat+(kk*NN)+kk);
		for (rows2=kk+1;rows2<NN;rows2++) { 
                       	    *(mat+(rows*NN)+rows2) = *(mat+(rows*NN)+rows2) - 
                       	    *(mat+(rows*NN)+kk) * *(mat+(kk*NN)+rows2) ;
                }
	}
	free(thread_args);
	pthread_exit(NULL);

}
