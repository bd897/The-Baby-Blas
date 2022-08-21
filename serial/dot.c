#ifdef __cplusplus
extern "C" {
#endif
	double dot_( int *threads, int *len, double *va, double *vb);
#ifdef __cplusplus
}
#endif

/* Serial Dot */

// Computes the tensor product of two vectors
double dot_( int *threads, int *len, double *va, double *vb){

	int i;
	int N = *len;
	double val = 0.0;

#ifdef OPT

	int mod;
	const int stride = 8;
	mod = N % stride;

	for(i=0; i<mod; i++){
		val += *(va+i) * *(vb+i);
	}
	for(i=mod; i<N; i+=stride){
		val +=   (*(va+i) * *(vb+i))
			+(*(va+i+1) * *(vb+i+1))
			+(*(va+i+2) * *(vb+i+2))
			+(*(va+i+3) * *(vb+i+3))
			+(*(va+i+4) * *(vb+i+4))
			+(*(va+i+5) * *(vb+i+5))
			+(*(va+i+6) * *(vb+i+6))
			+(*(va+i+7) * *(vb+i+7));
	}		

#else
	for(i=0; i<N; i++){
		val += (*(va+i) * *(vb+i));
	}
#endif
	return val;
}                                                                
