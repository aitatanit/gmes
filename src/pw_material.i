%module pw_material

%{
#define SWIG_FILE_WITH_INIT
#include "pw_material.hh"
#include "pw_dummy.hh"
#include "pw_const.hh"
#include "pw_dielectric.hh"
#include "pw_upml.hh"
#include "pw_cpml.hh"
#include "pw_drude.hh"
%}

%include "complex.i"
%include "numpy.i"
%numpy_typemaps(std::complex<double>, NPY_CDOUBLE, int)

%init %{
import_array();
%}

// Declare numpy typemaps.
%define %apply_numpy_typemaps(TYPE)
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const in_field1, int in1_dim1, int in1_dim2, int in1_dim3)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const in_field2, int in2_dim1, int in2_dim2, int in2_dim3)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3)};

%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const ex, int ex_x_size, int ex_y_size, int ex_z_size)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const ey, int ey_x_size, int ey_y_size, int ey_z_size)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const ez, int ez_x_size, int ez_y_size, int ez_z_size)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const hx, int hx_x_size, int hx_y_size, int hx_z_size)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const hy, int hy_x_size, int hy_y_size, int hy_z_size)};
%apply (TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(const TYPE* const hz, int hz_x_size, int hz_y_size, int hz_z_size)};
      
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const ex, int ex_x_size, int ex_y_size, int ex_z_size)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const ey, int ey_x_size, int ey_y_size, int ey_z_size)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const ez, int ez_x_size, int ez_y_size, int ez_z_size)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const hx, int hx_x_size, int hx_y_size, int hx_z_size)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const hy, int hy_x_size, int hy_y_size, int hy_z_size)};
%apply (TYPE* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
      {(TYPE* const hz, int hz_x_size, int hz_y_size, int hz_z_size)};
%enddef    /* apply_numpy_typemaps() macro */

%apply_numpy_typemaps(double)
%apply_numpy_typemaps(std::complex<double>)

%apply (int* IN_ARRAY1, int DIM1) {(const int* const idx, int size)};
%apply (double* IN_ARRAY1, int DIM1) {(const double * const omega_p, int omega_p_size)};
%apply (double* IN_ARRAY1, int DIM1) {(const double * const gamma_p, int gamma_p_size)};

// Declare the Pythonic interfaces.
%define %property(py, cpp, prop, get, set)
%feature("shadow") cpp::set %{ %}
%feature("shadow") cpp::get %{
__swig_setmethods__["prop"] = eval("_"+__name__.split('.')[-1]).##py##_##set
__swig_getmethods__["prop"] = eval("_"+__name__.split('.')[-1]).##py##_##get
if _newclass:prop = property(eval("_"+__name__.split('.')[-1]).##py##_##get, eval("_"+__name__.split('.')[-1]).##py##_##set)
%}
%enddef    /* property() macro */

%property(PwMaterialReal, gmes::PwMaterial<double>, i, get_i, set_i)
%property(PwMaterialReal, gmes::PwMaterial<double>, j, get_j, set_j)
%property(PwMaterialReal, gmes::PwMaterial<double>, k, get_k, set_k)

%property(PwMaterialCmplx, gmes::PwMaterial<std::complex<double> >, i, get_i, set_i)
%property(PwMaterialCmplx, gmes::PwMaterial<std::complex<double> >, j, get_j, set_j)
%property(PwMaterialCmplx, gmes::PwMaterial<std::complex<double> >, k, get_k, set_k)

%property(MaterialElectricReal, gmes::MaterialElectric<double>, epsilon, get_epsilon, set_epsilon)
%property(MaterialMagneticReal, gmes::MaterialMagnetic<double>, mu, get_mu, set_mu)

%property(MaterialElectricCmplx, gmes::MaterialElectric<std::complex<double> >, epsilon, get_epsilon, set_epsilon)
%property(MaterialMagneticCmplx, gmes::MaterialMagnetic<std::complex<double> >, mu, get_mu, set_mu)

%property(ConstElectricReal, gmes::ConstElectric<double>, value, get_value, set_value)
%property(ConstMagneticReal, gmes::ConstMagnetic<double>, value, get_value, set_value)

%property(ConstElectricCmplx, gmes::ConstElectric<std::complex<double> >, value, get_value, set_value)
%property(ConstMagneticCmplx, gmes::ConstMagnetic<std::complex<double> >, value, get_value, set_value)

// Include the header file to be wrapped
%include "pw_material.hh"
%include "pw_dummy.hh"
%include "pw_const.hh"
%include "pw_dielectric.hh"
%include "pw_upml.hh"
%include "pw_cpml.hh"
%include "pw_drude.hh"

// Instantiate template classes
%define %template_wrap(T, postfix) 
%template(PwMaterial ## postfix) gmes::PwMaterial<T >;
%template(MaterialElectric ## postfix) gmes::MaterialElectric<T >;
%template(MaterialMagnetic ## postfix) gmes::MaterialMagnetic<T >;

%template(DummyElectric ## postfix) gmes::DummyElectric<T >;
%template(DummyMagnetic ## postfix) gmes::DummyMagnetic<T >;
%template(DummyEx ## postfix) gmes::DummyEx<T >;
%template(DummyEy ## postfix) gmes::DummyEy<T >;
%template(DummyEz ## postfix) gmes::DummyEz<T >;
%template(DummyHx ## postfix) gmes::DummyHx<T >;
%template(DummyHy ## postfix) gmes::DummyHy<T >;
%template(DummyHz ## postfix) gmes::DummyHz<T >;

%template(ConstElectric ## postfix) gmes::ConstElectric<T >;
%template(ConstMagnetic ## postfix) gmes::ConstMagnetic<T >;
%template(ConstEx ## postfix) gmes::ConstEx<T >;
%template(ConstEy ## postfix) gmes::ConstEy<T >;
%template(ConstEz ## postfix) gmes::ConstEz<T >;
%template(ConstHx ## postfix) gmes::ConstHx<T >;
%template(ConstHy ## postfix) gmes::ConstHy<T >;
%template(ConstHz ## postfix) gmes::ConstHz<T >;

%template(DielectricElectric ## postfix) gmes::DielectricElectric<T >;
%template(DielectricMagnetic ## postfix) gmes::DielectricMagnetic<T >;
%template(DielectricEx ## postfix) gmes::DielectricEx<T >;
%template(DielectricEy ## postfix) gmes::DielectricEy<T >;
%template(DielectricEz ## postfix) gmes::DielectricEz<T >;
%template(DielectricHx ## postfix) gmes::DielectricHx<T >;
%template(DielectricHy ## postfix) gmes::DielectricHy<T >;
%template(DielectricHz ## postfix) gmes::DielectricHz<T >;

%template(UpmlElectric ## postfix) gmes::UpmlElectric<T >;
%template(UpmlMagnetic ## postfix) gmes::UpmlMagnetic<T >;
%template(UpmlEx ## postfix) gmes::UpmlEx<T >;
%template(UpmlEy ## postfix) gmes::UpmlEy<T >;
%template(UpmlEz ## postfix) gmes::UpmlEz<T >;
%template(UpmlHx ## postfix) gmes::UpmlHx<T >;
%template(UpmlHy ## postfix) gmes::UpmlHy<T >;
%template(UpmlHz ## postfix) gmes::UpmlHz<T >;

%template(CpmlElectric ## postfix) gmes::CpmlElectric<T >;
%template(CpmlMagnetic ## postfix) gmes::CpmlMagnetic<T >;
%template(CpmlEx ## postfix) gmes::CpmlEx<T >;
%template(CpmlEy ## postfix) gmes::CpmlEy<T >;
%template(CpmlEz ## postfix) gmes::CpmlEz<T >;
%template(CpmlHx ## postfix) gmes::CpmlHx<T >;
%template(CpmlHy ## postfix) gmes::CpmlHy<T >;
%template(CpmlHz ## postfix) gmes::CpmlHz<T >;

%template(DrudeElectric ## postfix) gmes::DrudeElectric<T >;
%template(DrudeEx ## postfix) gmes::DrudeEx<T >;
%template(DrudeEy ## postfix) gmes::DrudeEy<T >;
%template(DrudeEz ## postfix) gmes::DrudeEz<T >;
%template(DrudeHx ## postfix) gmes::DrudeHx<T >;
%template(DrudeHy ## postfix) gmes::DrudeHy<T >;
%template(DrudeHz ## postfix) gmes::DrudeHz<T >;
%enddef    /* template_wrap() macro */

%template_wrap(double, Real)
%template_wrap(std::complex<double>, Cmplx)
