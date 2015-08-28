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
#include "pw_lorentz.hh"
#include "pw_dcp.hh"
#include "pw_dm2.hh"
%}

%include <std_string.i>
%include <std_complex.i>
%include "numpy.i"

%numpy_typemaps(std::complex<double>, NPY_CDOUBLE, int)
%apply size_t { gmes::IdxCnt::size_type }; 

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

%apply (int* IN_ARRAY1, int DIM1) {(const int* const idx, int idx_size)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* const a, int a_size1, int a_size2)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* const b, int b_size1, int b_size2)};
%apply (std::complex<double>* IN_ARRAY2, int DIM1, int DIM2) {(const std::complex<double>* const b, int b_size1, int b_size2)};
%apply (double* IN_ARRAY1, int DIM1) {(const double* const c, int c_size)};

%apply (double* IN_ARRAY1, int DIM1) {(const double* const omega, int omega_size)};
%apply (double* IN_ARRAY1, int DIM1) {(const double* const n, int n_size)};

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* const u, int u_size)};
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* const v, int v_size)};
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* const w, int w_size)};

// Include the header file to be wrapped
%include "pw_material.hh"
%include "pw_dummy.hh"
%include "pw_const.hh"
%include "pw_dielectric.hh"
%include "pw_upml.hh"
%include "pw_cpml.hh"
%include "pw_drude.hh"
%include "pw_lorentz.hh"
%include "pw_dcp.hh"
%include "pw_dm2.hh"

// Instantiate template classes
%define %linear_wrap(T, postfix)

%template(ElectricParam ## postfix) gmes::ElectricParam<T >;
%template(MagneticParam ## postfix) gmes::MagneticParam<T >;
%template(PwMaterial ## postfix) gmes::PwMaterial<T >;
%template(MaterialElectric ## postfix) gmes::MaterialElectric<T >;
%template(MaterialMagnetic ## postfix) gmes::MaterialMagnetic<T >;

%template(DummyElectricParam ## postfix) gmes::DummyElectricParam<T >;
%template(DummyMagneticParam ## postfix) gmes::DummyMagneticParam<T >;
%template(DummyElectric ## postfix) gmes::DummyElectric<T >;
%template(DummyMagnetic ## postfix) gmes::DummyMagnetic<T >;
%template(DummyEx ## postfix) gmes::DummyEx<T >;
%template(DummyEy ## postfix) gmes::DummyEy<T >;
%template(DummyEz ## postfix) gmes::DummyEz<T >;
%template(DummyHx ## postfix) gmes::DummyHx<T >;
%template(DummyHy ## postfix) gmes::DummyHy<T >;
%template(DummyHz ## postfix) gmes::DummyHz<T >;

%template(ConstElectricParam ## postfix) gmes::ConstElectricParam<T >;
%template(ConstMagneticParam ## postfix) gmes::ConstMagneticParam<T >;
%template(ConstElectric ## postfix) gmes::ConstElectric<T >;
%template(ConstMagnetic ## postfix) gmes::ConstMagnetic<T >;
%template(ConstEx ## postfix) gmes::ConstEx<T >;
%template(ConstEy ## postfix) gmes::ConstEy<T >;
%template(ConstEz ## postfix) gmes::ConstEz<T >;
%template(ConstHx ## postfix) gmes::ConstHx<T >;
%template(ConstHy ## postfix) gmes::ConstHy<T >;
%template(ConstHz ## postfix) gmes::ConstHz<T >;

// Non-dispersive linear isotropic dielectrics
%template(DielectricElectricParam ## postfix) gmes::DielectricElectricParam<T >;
%template(DielectricMagneticParam ## postfix) gmes::DielectricMagneticParam<T >;
%template(DielectricElectric ## postfix) gmes::DielectricElectric<T >;
%template(DielectricMagnetic ## postfix) gmes::DielectricMagnetic<T >;
%template(DielectricEx ## postfix) gmes::DielectricEx<T >;
%template(DielectricEy ## postfix) gmes::DielectricEy<T >;
%template(DielectricEz ## postfix) gmes::DielectricEz<T >;
%template(DielectricHx ## postfix) gmes::DielectricHx<T >;
%template(DielectricHy ## postfix) gmes::DielectricHy<T >;
%template(DielectricHz ## postfix) gmes::DielectricHz<T >;

// UPML
%template(UpmlElectricParam ## postfix) gmes::UpmlElectricParam<T >;
%template(UpmlMagneticParam ## postfix) gmes::UpmlMagneticParam<T >;
%template(UpmlElectric ## postfix) gmes::UpmlElectric<T >;
%template(UpmlMagnetic ## postfix) gmes::UpmlMagnetic<T >;
%template(UpmlEx ## postfix) gmes::UpmlEx<T >;
%template(UpmlEy ## postfix) gmes::UpmlEy<T >;
%template(UpmlEz ## postfix) gmes::UpmlEz<T >;
%template(UpmlHx ## postfix) gmes::UpmlHx<T >;
%template(UpmlHy ## postfix) gmes::UpmlHy<T >;
%template(UpmlHz ## postfix) gmes::UpmlHz<T >;

// CPML
%template(CpmlElectricParam ## postfix) gmes::CpmlElectricParam<T >;
%template(CpmlMagneticParam ## postfix) gmes::CpmlMagneticParam<T >;
%template(CpmlElectric ## postfix) gmes::CpmlElectric<T >;
%template(CpmlMagnetic ## postfix) gmes::CpmlMagnetic<T >;
%template(CpmlEx ## postfix) gmes::CpmlEx<T >;
%template(CpmlEy ## postfix) gmes::CpmlEy<T >;
%template(CpmlEz ## postfix) gmes::CpmlEz<T >;
%template(CpmlHx ## postfix) gmes::CpmlHx<T >;
%template(CpmlHy ## postfix) gmes::CpmlHy<T >;
%template(CpmlHz ## postfix) gmes::CpmlHz<T >;

// Drude model
%template(DrudeElectricParam ## postfix) gmes::DrudeElectricParam<T >;
%template(DrudeMagneticParam ## postfix) gmes::DrudeMagneticParam<T >;
%template(DrudeElectric ## postfix) gmes::DrudeElectric<T >;
%template(DrudeEx ## postfix) gmes::DrudeEx<T >;
%template(DrudeEy ## postfix) gmes::DrudeEy<T >;
%template(DrudeEz ## postfix) gmes::DrudeEz<T >;
%template(DrudeHx ## postfix) gmes::DrudeHx<T >;
%template(DrudeHy ## postfix) gmes::DrudeHy<T >;
%template(DrudeHz ## postfix) gmes::DrudeHz<T >;

%extend gmes::DrudeElectricParam<T >
{
  void set(const double* const a, int a_size1, int a_size2,
	   const double* const c, int c_size)
  {
    for (int i = 0; i < a_size1; i++) {
      std::array<double, 3> tmp;
      std::copy(a + i * a_size2, a + i * a_size2 + 3, tmp.begin());
      $self->a.push_back(tmp);
    }
    
    std::copy(c, c + c_size, $self->c.begin());
    
    $self->q_now.resize(a_size1, T(0));
    $self->q_new.resize(a_size1, T(0));
  }
};

// Lerentz model
%template(LorentzElectricParam ## postfix) gmes::LorentzElectricParam<T >;
%template(LorentzMagneticParam ## postfix) gmes::LorentzMagneticParam<T >;
%template(LorentzElectric ## postfix) gmes::LorentzElectric<T >;
%template(LorentzEx ## postfix) gmes::LorentzEx<T >;
%template(LorentzEy ## postfix) gmes::LorentzEy<T >;
%template(LorentzEz ## postfix) gmes::LorentzEz<T >;
%template(LorentzHx ## postfix) gmes::LorentzHx<T >;
%template(LorentzHy ## postfix) gmes::LorentzHy<T >;
%template(LorentzHz ## postfix) gmes::LorentzHz<T >;

%extend gmes::LorentzElectricParam<T >
{
  void set(const double* const a, int a_size1, int a_size2,
	   const double* const c, int c_size)
  {
    for (int i = 0; i < a_size1; i++) {
      std::array<double, 3> tmp;
      std::copy(a + i * a_size2, a + i * a_size2 + 3, tmp.begin());
      $self->a.push_back(tmp);
    }

    std::copy(c, c + c_size, $self->c.begin());
    
    $self->l_now.resize(a_size1);
    $self->l_new.resize(a_size1);
  }
};

// ADE implementation of the Drude-critical points model
%template(DcpAdeElectricParam ## postfix) gmes::DcpAdeElectricParam<T >;
%template(DcpAdeMagneticParam ## postfix) gmes::DcpAdeMagneticParam<T >;
%template(DcpAdeElectric ## postfix) gmes::DcpAdeElectric<T >;
%template(DcpAdeEx ## postfix) gmes::DcpAdeEx<T >;
%template(DcpAdeEy ## postfix) gmes::DcpAdeEy<T >;
%template(DcpAdeEz ## postfix) gmes::DcpAdeEz<T >;
%template(DcpAdeHx ## postfix) gmes::DcpAdeHx<T >;
%template(DcpAdeHy ## postfix) gmes::DcpAdeHy<T >;
%template(DcpAdeHz ## postfix) gmes::DcpAdeHz<T >;

%extend gmes::DcpAdeElectricParam<T >
{
  void set(const double* const a, int a_size1, int a_size2,
	   const double* const b, int b_size1, int b_size2,
	   const double* const c, int c_size)
  {
    for (int i = 0; i < a_size1; i++) {
      std::array<double, 3> tmp;
      std::copy(a + i * a_size2, a + i * a_size2 + 3, tmp.begin());
      $self->a.push_back(tmp);
    }

    for (int i = 0; i < b_size1; i++) {
      std::array<double, 5> tmp;
      std::copy(b + i * b_size2, b + i * b_size2 + 5, tmp.begin());
      $self->b.push_back(tmp);
    }

    std::copy(c, c + c_size, $self->c.begin());
    
    $self->q_old.resize(a_size1);
    $self->q_now.resize(a_size1);
    $self->p_old.resize(b_size1);
    $self->p_now.resize(b_size1);
  }
};

// PLRC implementation of the Drude-critical points model
%template(DcpPlrcElectricParam ## postfix) gmes::DcpPlrcElectricParam<T >;
%template(DcpPlrcMagneticParam ## postfix) gmes::DcpPlrcMagneticParam<T >;
%template(DcpPlrcElectric ## postfix) gmes::DcpPlrcElectric<T >;
%template(DcpPlrcEx ## postfix) gmes::DcpPlrcEx<T >;
%template(DcpPlrcEy ## postfix) gmes::DcpPlrcEy<T >;
%template(DcpPlrcEz ## postfix) gmes::DcpPlrcEz<T >;
%template(DcpPlrcHx ## postfix) gmes::DcpPlrcHx<T >;
%template(DcpPlrcHy ## postfix) gmes::DcpPlrcHy<T >;
%template(DcpPlrcHz ## postfix) gmes::DcpPlrcHz<T >;

%extend gmes::DcpPlrcElectricParam<T >
{
  void set(const double* const a, int a_size1, int a_size2,
	   const std::complex<double>* const b, int b_size1, int b_size2,
	   const double* const c, int c_size)
  {
    for (int i = 0; i < a_size1; i++) {
      std::array<double, 3> tmp;
      std::copy(a + i * a_size2, a + i * a_size2 + 3, tmp.begin());
      $self->a.push_back(tmp);
    }

    for (int i = 0; i < b_size1; i++) {
      std::array<std::complex<double>, 3> tmp;
      std::copy(b + i * b_size2, b + i * b_size2 + 3, tmp.begin());
      $self->b.push_back(tmp);
    }

    std::copy(c, c + c_size, $self->c.begin());
    
    $self->psi_dp_re.resize(a_size1);
    $self->psi_dp_im.resize(a_size1);
    $self->psi_cp_re.resize(b_size1);
    $self->psi_cp_im.resize(b_size1);
  }
};

%enddef    /* linear_wrap() macro */

%define %nonlinear_wrap(T, postfix)

// density matrix implementation
%template(Dm2ElectricParam ## postfix) gmes::Dm2ElectricParam<T >;
%template(Dm2MagneticParam ## postfix) gmes::Dm2MagneticParam<T >;
%template(Dm2Electric ## postfix) gmes::Dm2Electric<T >;
%template(Dm2Ex ## postfix) gmes::Dm2Ex<T >;
%template(Dm2Ey ## postfix) gmes::Dm2Ey<T >;
%template(Dm2Ez ## postfix) gmes::Dm2Ez<T >;
%template(Dm2Hx ## postfix) gmes::Dm2Hx<T >;
%template(Dm2Hy ## postfix) gmes::Dm2Hy<T >;
%template(Dm2Hz ## postfix) gmes::Dm2Hz<T >;

%extend gmes::Dm2ElectricParam<T >
{
  void set(const double* const omega, int omega_size,
           const double* const n, int n_size)
  {
    for (int i = 0; i < omega_size; i++) {
      $self->omega.push_back(*(omega + i));
      $self->n_atom.push_back(*(n + i));
      
      std::array<double, 3> u_tmp;
      u_tmp.fill(static_cast<T>(0));
      $self->u.push_back(u_tmp);
    }
  }
};

%enddef    /* dm2_wrap() macro */

%linear_wrap(double, Real)
%linear_wrap(std::complex<double>, Cmplx)

%nonlinear_wrap(double, Real)
