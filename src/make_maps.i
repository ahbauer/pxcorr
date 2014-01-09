%module make_maps

%include exception.i
%include "stl.i"

%{
#define SWIG_FILE_WITH_INIT
extern void make_maps( const char *catalog_filename_c, const char *mask_filename_c, float *ang_means_c, int n_ang_bins0, float *ang_widths_c, int n_ang_bins, bool use_counts, bool use_mags, char *suffix_c );
%}

%include "numpy.i"
%init %{
import_array();
%}

%exception {
    try {
        $action
    } 
    catch( const std::exception & e ) {
        SWIG_exception(SWIG_RuntimeError, (std::string("C++ exception: ") + e.what()).c_str());
    }
}

%apply (float* IN_ARRAY1, int DIM1) {(float *ang_means_c, int n_ang_bins0)};
%apply (float* IN_ARRAY1, int DIM1) {(float *ang_widths_c, int n_ang_bins)};
extern void make_maps( const char *catalog_filename_c, const char *mask_filename_c, float *ang_means_c, int n_ang_bins0, float *ang_widths_c, int n_ang_bins, bool use_counts, bool use_mags, char *suffix_c );
