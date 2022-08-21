
#ifdef __cplusplus
extern "C" {
#endif
    void mvv_( int *threads, int *len, double *mat, double *vec, double *rvec);
#ifdef __cplusplus
    }
#endif

/* OPENMP MVV */

#include <omp.h>

// Computes the tensor product of two vectors 
void  mvv_( int *threads, int *len, double *mat, double *vec, double *rvec){

	int i, j;
	int N = *len;
	// Vars for loop unrolling
	int mod;
	int stride = 8;
	mod = N%stride;
	
	// Set the number of threads
 	omp_set_num_threads(*threads);

    #pragma omp parallel shared( i, N, rvec, mat, mod, stride)
    {
	#pragma omp for schedule(static)
	for(i=0; i<N; i++) *(rvec+i) = 0.0;	

	#pragma omp for private(j) 
	for (i=0; i<N; i++){
		for(j=0; j<mod; j++){
			*(rvec+i) += ( *(mat+((N*i)+j)) * *(vec+j));
		}
		for (j=mod; j<N; j+=stride){
			*(rvec+i) += ( ( *(mat+((N*i)+j)) * *(vec+j) )
				+ ( *(mat+((N*i)+j+1)) * *(vec+j+1))
				+ ( *(mat+((N*i)+j+2)) * *(vec+j+2))
				+ ( *(mat+((N*i)+j+3)) * *(vec+j+3))
				+ ( *(mat+((N*i)+j+4)) * *(vec+j+4))
				+ ( *(mat+((N*i)+j+5)) * *(vec+j+5))
				+ ( *(mat+((N*i)+j+6)) * *(vec+j+6))
				+ ( *(mat+((N*i)+j+7)) * *(vec+j+7)) );
        	}
	}
    }
	// Our resultant vector should be in rvec
}
