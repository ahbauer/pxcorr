%module correlate

%include exception.i
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
extern void correlate( char* mapn1, char* mapn2, int r, char* outfilename );
%}

%exception {
    try {
        $action
    } 
    catch( const std::exception & e ) {
        SWIG_exception(SWIG_RuntimeError, (std::string("C++ exception: ") + e.what()).c_str());
    }
    catch( ... ) {
        SWIG_exception(SWIG_RuntimeError, (std::string("C++ non-standard exception!")).c_str());
    }
}

%include "numpy.i"
%init %{
import_array();
%}

extern void correlate( char* mapn1, char* mapn2, int r, char* outfilename );
