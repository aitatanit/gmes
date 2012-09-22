%module constant

%{
#define SWIG_FILE_WITH_INIT
#include "constant.hh"
%}

// Get the NumPy typemaps
%include "numpy.i"

%init %{
import_array();
%}

%apply (double ARGOUT_ARRAY1[ANY]) {(double vector[3])};

%rename(_component) gmes::Component;
%rename(_electric) gmes::Electric;
%rename(_ex) gmes::Ex;
%rename(_ey) gmes::Ey;
%rename(_ez) gmes::Ez;
%rename(_magnetic) gmes::Magnetic;
%rename(_hx) gmes::Hx;
%rename(_hy) gmes::Hy;
%rename(_hz) gmes::Hz;
%rename(_electriccurrent) gmes::ElectricCurrent;
%rename(_jx) gmes::Jx;
%rename(_jy) gmes::Jy;
%rename(_jz) gmes::Jz;
%rename(_magneticcurrent) gmes::MagneticCurrent;
%rename(_mx) gmes::Mx;
%rename(_my) gmes::My;
%rename(_mz) gmes::Mz;
%rename(_directional) gmes::Directional;
%rename(_x) gmes::X;
%rename(_y) gmes::Y;
%rename(_z) gmes::Z;
%rename(_plus_x) gmes::PlusX;
%rename(_plus_y) gmes::PlusY;
%rename(_plus_z) gmes::PlusZ;
%rename(_minus_x) gmes::MinusX;
%rename(_minus_y) gmes::MinusY;
%rename(_minus_z) gmes::MinusZ;

// Include the header file to be wrapped
%include "constant.hh"

%pythoncode %{
class Component(_component):
    tag = _component_get_tag()
	
class Electric(Component):
  tag = _electric_get_tag()
    
class Ex(Electric):
    tag = _ex_get_tag()

    @classmethod
    def str(cls):
        return 'ex'

class Ey(Electric):
    tag = _ey_get_tag()

    @classmethod
    def str(cls):
        return 'ey'

class Ez(Electric):
    tag = _ez_get_tag()

    @classmethod
    def str(cls):
        return 'ez'

class Magnetic(Component):
    tag = _magnetic_get_tag()
    
class Hx(Magnetic):
    tag = _hx_get_tag()

    @classmethod
    def str(cls):
        return 'hx'

class Hy(Magnetic):
    tag = _hy_get_tag()

    @classmethod
    def str(cls):
        return 'hy'

class Hz(Magnetic):
    tag = _hz_get_tag()

    @classmethod
    def str(cls):
        return 'hz'

class ElectricCurrent(Component):
  tag = _electriccurrent_get_tag()
    
class Jx(ElectricCurrent):
    tag = _jx_get_tag()

    @classmethod
    def str(cls):
        return 'jx'

class Jy(ElectricCurrent):
    tag = _jy_get_tag()

    @classmethod
    def str(cls):
        return 'jy'

class Jz(ElectricCurrent):
    tag = _jz_get_tag()

    @classmethod
    def str(cls):
        return 'jz'

class MagneticCurrent(Component):
    tag = _magneticcurrent_get_tag()
    
class Mx(Magnetic):
    tag = _mx_get_tag()

    @classmethod
    def str(cls):
        return 'mx'

class My(MagneticCurrent):
    tag = _my_get_tag()

    @classmethod
    def str(cls):
        return 'my'

class Mz(MagneticCurrent):
    tag = _mz_get_tag()

    @classmethod
    def str(cls):
        return 'mz'

class Directional(_directional):
    tag = _directional_get_tag()
    
class X(Directional):
    tag = _x_get_tag()

    @classmethod
    def str(cls):
        return 'x'

class Y(Directional):
    tag = _y_get_tag()

    @classmethod
    def str(cls):
        return 'y'
    
class Z(Directional):
    tag = _z_get_tag()

    @classmethod
    def str(cls):
        return 'z'
    
class PlusX(X):
    tag = _plus_x_get_tag()
    vector = _plus_x_get_vector()
    @classmethod
    def str(cls):
        return '+x'

class MinusX(X):
    tag = _minus_x_get_tag()
    vector = _minus_x_get_vector()

    @classmethod
    def str(cls):
        return '-x'
    
class PlusY(Y):
    tag = _plus_y_get_tag()
    vector = _plus_y_get_vector()

    @classmethod
    def str(cls):
        return '+y'
        
class MinusY(Y):
    tag = _minus_y_get_tag()
    vector = _minus_y_get_vector()

    @classmethod
    def str(cls):
        return '-y'
    
class PlusZ(Z):
    tag = _plus_z_get_tag()
    vector = _plus_z_get_vector()

    @classmethod
    def str(cls):
        return '+z'
        
class MinusZ(Z):
    tag = _minus_z_get_tag()
    vector = _minus_z_get_vector()

    @classmethod
    def str(cls):
        return '-z'
%}
