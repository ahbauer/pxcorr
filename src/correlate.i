%module correlate
%{
#define SWIG_FILE_WITH_INIT
extern void correlate( char* mapn1, char* mapn2, char* sfx, int r, double **outarray, int *nout );
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double **outarray, int* nout)};
extern void correlate( char* mapn1, char* mapn2, char* sfx, int r, double **outarray, int *nout );
