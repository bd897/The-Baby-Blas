#ifdef __cplusplus
extern "C" {
#endif
    void vvm_(  int *threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
    }
#endif

/* Serial VVM */

// Computes the tensor product of two vectors 
void  vvm_( int *threads, int *len, double *va, double *vb, double *ma){

	int i, j;
	int N = *len;

// Unrolls Loop and uses temp
#ifdef OPT

	double tmp;
	int mod;
	const int stride = 8;
	mod = N%stride;

	for(i=0; i<N; i++){
		tmp = *(va+i);
		for(j=0; j<mod; j++){ 
			*(ma+(i*N)+j) = tmp * *(vb+j)
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

#else

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			*(ma+(i*N)+j) = *(va+i) * *(vb+j);
		}
	}

#endif
}
