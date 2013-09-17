%module correlate3pt
%{
#define SWIG_FILE_WITH_INIT
extern void correlate3pt( char* mapn1, char* mapn2, char* mapn3, char* sfx );
%}

%include "numpy.i"
%init %{
import_array();
%}

extern void correlate3pt( char* mapn1, char* mapn2, char* mapn3, char* sfx );
