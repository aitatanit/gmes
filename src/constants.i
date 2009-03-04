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

%pythoncode %{
class Component(Component):
    tag = Component_get_tag()
    
class Electric(Electric):
    tag = Electric_get_tag()
    
class Ex(Ex):
    tag = Ex_get_tag()
    
class Ey(Ey):
    tag = Ey_get_tag()
    
class Ez(Ez):
    tag = Ez_get_tag()
    
class Magnetic(Magnetic):
    tag = Magnetic_get_tag()
    
class Hx(Hx):
    tag = Hx_get_tag()
    
class Hy(Hy):
    tag = Hy_get_tag()
    
class Hz(Hz):
    tag = Hz_get_tag()
    
#class Directional(Directional):
#    tag = Directional_get_tag()
    
class X(X):
    tag = X_get_tag()
    
class Y(Y):
    tag = Y_get_tag()
    
class Z(Z):
    tag = Z_get_tag()
    
class PlusX(PlusX):
    tag = PlusX_get_tag()
    vector = PlusX_get_vector()
    
class MinusX(MinusX):
    tag = MinusX_get_tag()
    vector = MinusX_get_vector()
    
class PlusY(PlusY):
    tag = PlusY_get_tag()
    vector = PlusY_get_vector()
    
class MinusY(MinusY):
    tag = MinusY_get_tag()
    vector = MinusY_get_vector()
    
class PlusZ(PlusZ):
    tag = PlusZ_get_tag()
    vector = PlusZ_get_vector()
    
class MinusZ(MinusZ):
    tag = MinusZ_get_tag()
    vector = MinusZ_get_vector()
%}
