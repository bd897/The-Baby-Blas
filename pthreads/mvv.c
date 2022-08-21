#ifdef __cplusplus
extern "C" {
#endif
	double mvv_( int *threads, int *len, double *ma, double *vec, double *vr);
#ifdef __cplusplus
}
#endif

#include <pthread.h>
#include <stdlib.h>

// Prototype
void *mvv_thread_worker();

// Struct to pass thread
struct args{
	int N;
	int startRow;
	int stopRow;
	double *Aptr;
	double *Bptr;
	double *Cptr;
};

// Computes the vector product between a matrix and vector
double mvv_( int *threads, int *len, double *ma, double *vec, double *vr){

 	
	int i,j;
	int numThreads = *threads;
	int N = *len;
	int *numRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;

	// If threads are the same as the vector length, then do MVV in serial
	// If not do in parallel
	if( N < numThreads ){
		// MVV in serial
		for (i=0; i<N; i++){
			*(vr+i) = 0.0;
			for (j=0; j<N; j++){
				*(vr+i) += *(ma+(N*i)+j) * *(vec+j) ;
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
	
		// Makes the args and create each thread
		stopRow = 0;
		for(i=0; i<numThreads; i++){
			{
				startRow=stopRow;
				stopRow=startRow+*(numRows+i);
				thread_args = (struct args *) malloc(sizeof(struct args));
				thread_args -> N = N;
				thread_args -> startRow = startRow;
				thread_args -> stopRow = stopRow;
				thread_args -> Aptr = ma;
				thread_args -> Bptr = vec;
				thread_args -> Cptr = vr;
				pthread_create( thread_id+i, NULL, &mvv_thread_worker, thread_args );
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


void *mvv_thread_worker( struct args *thread_args){

	int i, j, mod, stride;
	int rowStart, rowStop, NN;
	double *mat, *vec, *res;

	// Grab data from struct
	NN       = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop  = thread_args->stopRow;
	mat      = thread_args->Aptr;
	vec      = thread_args->Bptr;
	res      = thread_args->Cptr;

	stride = 8;
	mod = NN%stride;

	// MVV Algorithm
	for (i=rowStart; i<rowStop; i++){
		 *(res+i) = 0.0;
		for(j=0; j<mod; j++){
			*(res+i) += ( *(mat+((NN*i)+j)) * *(vec+j));
		}
		for(j=mod; j<NN; j+=stride){
			*(res+i) = *(res+i) + ( ( *(mat+((NN*i)+j)) * *(vec+j))
		 		+ ( *(mat+((NN*i)+j+1)) * *(vec+j+1))
				+ ( *(mat+((NN*i)+j+2)) * *(vec+j+2))
				+ ( *(mat+((NN*i)+j+3)) * *(vec+j+3))
				+ ( *(mat+((NN*i)+j+4)) * *(vec+j+4))
				+ ( *(mat+((NN*i)+j+5)) * *(vec+j+5))
				+ ( *(mat+((NN*i)+j+6)) * *(vec+j+6))
				+ ( *(mat+((NN*i)+j+7)) * *(vec+j+7)) );
        	}
	}

	// Clean Up
	free(thread_args);	
	// Send the thred back
	pthread_exit(NULL);
}
