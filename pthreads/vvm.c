#ifdef __cplusplus
extern "C" {
#endif
	double vvm_( int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

// Prototype
void *vvm_thread_worker();

// Struct to pass thread
struct args{
	int N;
	int startRow;
	int stopRow;
	double *Aptr;
	double *Bptr;
	double *Cptr;
};

// Computes the tensor product of two vectors
double vvm_( int *threads, int *len, double *va, double *vb, double *ma){
 	
	int i, j, mod, stride, startRow, stopRow; 
	int numThreads = *threads;
	int N = *len;
	int *numRows;
	pthread_t *thread_id;
	struct args *thread_args;
	double tmp;

	// If threads are the same as the vector length, then do normal dot
	// If else, do in parallel
	if( N < numThreads ){
		stride = 8;
		mod = N%stride;	

		for(i=0; i<N; i++){
			tmp = *(va+i);
			for(j=0; j<mod; j++){ 
				*(ma+(i*N)+j) = tmp * *(vb+j);
			}
			for(j=mod; j<N; j+=stride){
				*(ma+(i*N)+j) = tmp * *(vb+j);
				*(ma+(i*N)+j+1) = tmp * *(vb+j+1);
				*(ma+(i*N)+j+2) = tmp * *(vb+j+2);
				*(ma+(i*N)+j+3) = tmp * *(vb+j+3);
				*(ma+(i*N)+j+4) = tmp * *(vb+j+4);
				*(ma+(i*N)+j+5) = tmp * *(vb+j+5);
				*(ma+(i*N)+j+6) = tmp * *(vb+j+6);
				*(ma+(i*N)+j+7) = tmp * *(vb+j+7);
	       		}
		}
	}
	else{
		// Allocate memory for thread id's and numRows
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
		numRows = (int *) malloc (numThreads * sizeof(int));

		// Divide up elements for thread to work
		for(i=0; i<numThreads; i++){
			*(numRows+i) = N/numThreads;
		}
		for(i=0; i< N%numThreads; i++){
			*(numRows+i) = *(numRows+i) + 1;
		}
	
		// Makes the args and creates each thread
		stopRow = 0;
		for(i=0; i<numThreads; i++){
			{
				startRow=stopRow;
				stopRow=startRow+*(numRows+i);
				thread_args = (struct args *) malloc(sizeof(struct args));
				thread_args -> N = N;
				thread_args -> startRow = startRow;
				thread_args -> stopRow = stopRow;
				thread_args -> Aptr = va;
				thread_args -> Bptr = vb;
				thread_args -> Cptr = ma;
				pthread_create( thread_id+i, NULL, &vvm_thread_worker, thread_args );
			}
		}

		// Retrieve info from threads
		for(i=0; i<numThreads; i++){
			pthread_join(*(thread_id+i), NULL);
		}
	
		// Clean Up
		free(numRows);
		free(thread_id);
	}
}


void *vvm_thread_worker( struct args *thread_args){

	int i, j, NN, mod, stride, rowStart, rowStop, N;
	double tmp;
	double *va, *vb, *ma;

	// Grab data from struct
	NN        = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop  = thread_args->stopRow;
	va       = thread_args->Aptr;
	vb       = thread_args->Bptr;
	ma       = thread_args->Cptr;

	stride = 8;
	mod = NN%stride;

	for(i=rowStart; i<rowStop; i++){
		tmp = *(va+i);
		for(j=0; j<mod; j++){ 
			*(ma+(i*NN)+j) = tmp * *(vb+j);
		}
		for(j=mod; j<NN; j+=stride){
			*(ma+(i*NN)+j) = tmp * *(vb+j);
			*(ma+(i*NN)+j+1) = tmp * *(vb+j+1);
			*(ma+(i*NN)+j+2) = tmp * *(vb+j+2);
			*(ma+(i*NN)+j+3) = tmp * *(vb+j+3);
			*(ma+(i*NN)+j+4) = tmp * *(vb+j+4);
			*(ma+(i*NN)+j+5) = tmp * *(vb+j+5);
			*(ma+(i*NN)+j+6) = tmp * *(vb+j+6);
			*(ma+(i*NN)+j+7) = tmp * *(vb+j+7);
	       	}
	}

	// Clean Up
	free(thread_args);
	// Send the thread back
	pthread_exit(NULL);
}





















