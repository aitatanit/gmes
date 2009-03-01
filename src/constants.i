%module constants

%{
#define SWIG_FILE_WITH_INIT
#include "constants.hh"
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
import_array();
%}

%apply (double ARGOUT_ARRAY1[ANY]) {(double vector[3])};

// Include the header file to be wrapped
%include "constants.hh"
