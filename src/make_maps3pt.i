%module make_maps3pt
%{
#define SWIG_FILE_WITH_INIT
extern void make_maps3pt( const char *catalog_filename_c, const char *mask_filename_c, float *angles_c, int n_angles, float *r12s_c, int n_r12s, float *r23s_c, int n_r23s, char *suffix_c );

%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (float* IN_ARRAY1, int DIM1) {(float *angles_c, int n_angles)};
%apply (float* IN_ARRAY1, int DIM1) {(float *r12s_c, int n_r12s)};
%apply (float* IN_ARRAY1, int DIM1) {(float *r23s_c, int n_r23s)};
extern void make_maps3pt( const char *catalog_filename_c, const char *mask_filename_c, float *angles_c, int n_angles, float *r12s_c, int n_r12s, float *r23s_c, int n_r23s, char *suffix_c );
