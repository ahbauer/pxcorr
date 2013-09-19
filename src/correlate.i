%module correlate
%{
#define SWIG_FILE_WITH_INIT
extern void correlate( char* mapn1, char* mapn2, char* sfx, int r );
%}

%include "numpy.i"
%init %{
import_array();
%}

extern void correlate( char* mapn1, char* mapn2, char* sfx, int r );
