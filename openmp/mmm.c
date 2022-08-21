#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
    }
#endif

/*  O p e n M P     C O D E  */

#include <omp.h>

void mmm_( int *threads, int *len,  double *a, double *b, double *c ){

    	int i, j, k;
    	int N = *len;

	// Set the number of threads to use here
	omp_set_num_threads(*threads);

	#pragma omp parallel for shared(N) private(i,j,k) schedule(static) collapse(2)
	for (i=0; i<N; i++) {
        	for (j=0; j<N; j++) {
            		*(c+(i*N)+j) = 0.0;
            		for (k=0;k<N;k++){
                		*(c+(i*N)+j) += *(a+(i*N)+k) * *(b+(k*N)+j); 
            		}
        	}
    	}
}
