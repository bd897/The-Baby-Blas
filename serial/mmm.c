#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double *c);
#ifdef __cplusplus
    }
#endif

/* Serial MMM */

void mmm_( int *threads, int *len,  double *a, double *b, double *c ){

	int i, j, k;
	int N = *len;

#ifdef STRIP
	int mod;
	const int stride = 8;
	mod = N % stride;

	for (i=0; i<N; i++) {
        	for (j=0; j<N; j++) {
            		*(c+(i*N+j)) = 0.0;
            		for (k=0;k<mod;k++){
                		*(c+(i*N+j)) += *(a+(i*N+k)) * *(b+(k*N+j)); 
            		}
            		for (k=mod;k<N;k+=stride) {
                		*(c+(i*N+j)) += *(a+(i*N)+k) * *(b+(k*N)+j) 
                                   	+ *(a+(i*N+k+1)) * *(b+((k+1)*N+j)) 
                                   	+ *(a+(i*N+k+2)) * *(b+((k+2)*N+j)) 
                                   	+ *(a+(i*N+k+3)) * *(b+((k+3)*N+j)) 
                                   	+ *(a+(i*N+k+4)) * *(b+((k+4)*N+j)) 
                                   	+ *(a+(i*N+k+5)) * *(b+((k+5)*N+j)) 
                                   	+ *(a+(i*N+k+6)) * *(b+((k+6)*N+j)) 
                                   	+ *(a+(i*N+k+7)) * *(b+((k+7)*N+j)); 
            		}
        	}
	}
#else
	for (i=0; i<N; i++) {
        	for (j=0; j<N; j++) {
            		*(c+(i*N)+j) = 0.0;
            		for (k=0;k<N;k++){
                		*(c+(i*N)+j) += *(a+(i*N)+k) * *(b+(k*N)+j); 
			}
		}	
	}
#endif
}


