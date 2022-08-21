
#ifdef __cplusplus
extern "C" {
#endif
    void mvv_( int *threads, int *len, double *ma, double *vec, double *rvec);
#ifdef __cplusplus
    }
#endif

/* Serial MVV */

// Computes the vector product between a matrix and a vector
void  mvv_( int *threads, int *len, double *mat, double *vec, double *rvec){

	int i, j;
	int N = *len;

// Unrolls Loop
#ifdef OPT

	int mod;
	const int stride = 8;
	mod = N%stride;

	for(i=0; i<N; i++){
		*(rvec+i) = 0.0;
		for(j=0; j<mod; j++){
			*(rvec+i) += ( *(mat+((N*i)+j)) * *(vec+j));
		}
		for(j=mod; j<N; j+=stride){
			*(rvec+i) = *(rvec+i) + ( ( *(mat+((N*i)+j)) * *(vec+j))
		 		+ ( *(mat+((N*i)+j+1)) * *(vec+j+1))
				+ ( *(mat+((N*i)+j+2)) * *(vec+j+2))
				+ ( *(mat+((N*i)+j+3)) * *(vec+j+3))
				+ ( *(mat+((N*i)+j+4)) * *(vec+j+4))
				+ ( *(mat+((N*i)+j+5)) * *(vec+j+5))
				+ ( *(mat+((N*i)+j+6)) * *(vec+j+6))
				+ ( *(mat+((N*i)+j+7)) * *(vec+j+7)) );
        	}
	}
#else
	for (i=0; i<N; i++){
		*(rvec+i) = 0.0;
		for (j=0; j<N; j++){
			*(rvec+i) = *(rvec+i) + ( *(mat+((N*i)+j)) * *(vec+j) );
		}
	}
#endif
}
