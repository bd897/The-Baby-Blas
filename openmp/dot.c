#ifdef __cplusplus
extern "C" {
#endif
	double dot_( int *threads, int *len, double *va, double *vb);
#ifdef __cplusplus
}
#endif

/* DOT IN OPENMP */

#include <omp.h>

// Computes the tensor product of two vectors
double dot_( int *threads, int *len, double *va, double *vb){

	int i;
	int N = *len;
	double ans = 0.0;
	
	// Sets the Threads
	omp_set_num_threads(*threads);	

	#pragma omp parallel for reduction(+:ans)
	for(i=0; i<N; i++){
		 ans += (*(va+i) * *(vb+i));
	}

	return ans;
}                                                                
