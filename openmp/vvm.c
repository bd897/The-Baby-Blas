#ifdef __cplusplus
extern "C" {
#endif
    void vvm_( int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
    }
#endif

/* OPENMP VVM */

#include <omp.h>

// Computes the tensor product of two vectors 
void  vvm_( int *threads, int *len, double *va, double *vb, double *ma){

	int i, j, stride, mod;
	int numThreads = *threads;
	int N = *len;
	double tmp;

	// Sets number of threads
	omp_set_num_threads(numThreads);

	stride = 8;
	mod = N%stride;	

	#pragma omp parallel for private(j,tmp) schedule(static)
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



