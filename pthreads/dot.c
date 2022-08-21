#ifdef __cplusplus
extern "C" {
#endif
	double dot_( int *threads, int *len, double *va, double *vb);
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

// Prototype
void *dot_thread_worker();

// Struct to pass thread
struct args{
	int N;
	int startRow;
	int stopRow;
	double *Aptr;
	double *Bptr;
	double *sum;
};

// Create a mutex variable for our value
pthread_mutex_t val_mutex;


// Computes the dot product of two vectors
double dot_( int *threads, int *len, double *va, double *vb){
 	
	int numThreads = *threads;
	int vlen = *len;
	int *numRows;
	int startRow, stopRow;
	pthread_t *thread_id;
	struct args *thread_args;
	double *sum;
	double val;

	// Initialize our total sum
	val = 0.0;
	
	// If the number of threads are the same as the vec length, then do normal dot
	// If not do in parallel
	if( vlen < numThreads ){
		for(int i=0; i<vlen; i++){
			val += (*(va+i) * *(vb+i));
		}
	}
	else{
		// Allocate memory to thread id, numRows, and our sum
		thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));
		numRows = (int *) malloc (numThreads * sizeof(int));
		sum = (double *) malloc (sizeof(double));

		// Initialize our sum
		*sum = 0.0;
		
		// Divide up elements for thread to work
		for(int i=0; i<numThreads; i++){
			*(numRows+i) = vlen/numThreads;
		}
		for(int i=0; i< vlen%numThreads; i++){
			*(numRows+i) = *(numRows+i) + 1;
		}

		// Makes the args and creates each thread
		stopRow = 0;
		for(int i=0; i<numThreads; i++){
			{
				startRow=stopRow;
				stopRow=startRow+*(numRows+i);
				thread_args = (struct args *) malloc(sizeof(struct args));
				thread_args -> N = vlen;
				thread_args -> startRow = startRow;
				thread_args -> stopRow = stopRow;
				thread_args -> Aptr = va;
				thread_args -> Bptr = vb;				
				thread_args -> sum = sum;				
				pthread_create(thread_id+i,NULL,&dot_thread_worker,thread_args);
			}
		}
		// Collect the threads
		for(int i=0; i<numThreads; i++){
			pthread_join(*(thread_id+i), NULL);
		}
		
		// Retrieve answer
		val = *sum;
		
		// Clean Up
		free(numRows);
		free(thread_id);
		free(sum);
	}
	// Return our dot product
	return val;
}


void *dot_thread_worker( struct args *thread_args){

	int i, N, rowStart, rowStop;
	double *va, *vb, *sum;	
	double partSum;	

	// Grab data from struct
	N        = thread_args->N;
	rowStart = thread_args->startRow;
	rowStop  = thread_args->stopRow;
	va       = thread_args->Aptr;
	vb       = thread_args->Bptr;
	sum      = thread_args->sum;	

	// Intialize the partial sum at 0
	partSum = 0.0;

	// Do the dot product
	for(i=rowStart; i<rowStop; i++){
		partSum = partSum + (*(va+i) * *(vb+i));
	}

	// Locks our total sum, to have the partial sum of thread added to total sum
	pthread_mutex_lock(&val_mutex);
	*sum += partSum;
	pthread_mutex_unlock(&val_mutex);	

	// Clean Up
	free(thread_args);	
	// Send the thread back
	pthread_exit(NULL);
}

