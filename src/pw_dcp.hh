/* The implementation of Drude-critical points model based on the 
 * following articles.
 *
 * J. Leng, J. Opsal, H. Chu, M. Senko, and D. E. Aspnes, 
 * "Analytic representations of the dielectric functions of 
 * materials for device and structural modeling," Thin Solid Films
 * 313-314, 132-136 (1998).
 *
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model
 * for the optical properties of gold," J. Chem. Phys. 125, 
 * 164705-3 (2006).
 */

#ifndef PW_DCP_HH_
#define PW_DCP_HH_

#include <array>
#include <complex>
#include <numeric>
#include <string>
#include <vector>
#include "pw_dielectric.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

// The classes should be rewritten using template specialization
// to increase the calculation speed.
namespace gmes
{
  /* The following auxiliary differential equation(ADE) 
   * implementation of the DCP model is based on the following 
   * article.
   *
   * L. J. Prokopeva, J. D. Borneman, and A. V. Kildishev, "Optical
   * Dispersion Models for Time-Domain Modeling of Metal-
   * Dielectric Nanostructures," IEEE Trans. Magn. 47, 1150-1153 
   * (2011).
   */

  typedef std::vector<std::array<double, 3> > AdeCoeffA;
  typedef std::vector<std::array<double, 5> > AdeCoeffB;
  typedef std::array<double, 4> AdeCoeffC;
    
  template <typename T> 
  struct DcpAdeElectricParam: ElectricParam<T>
  {
    // parameters for the ADE of the Drude terms
    AdeCoeffA a;
    // parameters for the ADE of critical point terms
    AdeCoeffB b;
    // parameters for the electric field update equations
    AdeCoeffC c;
    T e_old;
    std::vector<T> q_old, q_now, p_old, p_now;
  }; // template DcpAdeElectricParam
  
  template <typename T> 
  struct DcpAdeMagneticParam: MagneticParam<T>
  {
  }; // DcpAdemagneticParam

  template <typename T> 
  class DcpAdeElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpAdeElectric<T>::tag;
    }

    double
    get_eps_inf(const int* const idx, int idx_size) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      const int i = position(index);
      if (i < 0)
	return 0;
      else
	return param_list[i].eps_inf;
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      const auto& dcp_param = *static_cast<const DcpAdeElectricParam<T> * const>(pm_param_ptr);
      
      idx_list.push_back(index);
      param_list.push_back(dcp_param);

      return this;
    }
    
    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dcp_ptr = static_cast<const DcpAdeElectric<T>*>(pm_ptr);
      std::copy(dcp_ptr->idx_list.begin(), dcp_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(dcp_ptr->param_list.begin(), dcp_ptr->param_list.end(), std::back_inserter(param_list));
      return this;
    }

    T 
    dps_sum(const T& init, const DcpAdeElectricParam<T>& dcp_param) const
    {
      const auto& a = dcp_param.a;
      const auto& q_old = dcp_param.q_old;
      const auto& q_now = dcp_param.q_now;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	sum += (1 - a[i][1]) * q_now[i] - a[i][0] * q_old[i];
      }

      return sum;
    }
    
    T 
    cps_sum(const T& init, const DcpAdeElectricParam<T>& dcp_param) const
    {
      const auto& b = dcp_param.b;
      const auto& p_old = dcp_param.p_old;
      const auto& p_now = dcp_param.p_now;

      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	sum += (1 - b[i][1]) * p_now[i] - b[i][0] * p_old[i];
      }
      
      return sum;
    }

    void 
    update_q(const T& e_old, const T& e_now, const T& e_new,
	     DcpAdeElectricParam<T>& dcp_param)
    {
      const auto& a = dcp_param.a;
      auto& q_old = dcp_param.q_old;
      auto& q_now = dcp_param.q_now;

      std::vector<T> q_new(a.size());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * (e_old + 2.0 * e_now + e_new);
      }
      
      std::copy(q_now.begin(), q_now.end(), q_old.begin());
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
    }
    
    void 
    update_p(const T& e_old, const T& e_now, const T& e_new,
	     DcpAdeElectricParam<T>& dcp_param)
    {
      const auto& b = dcp_param.b;
      auto& p_old = dcp_param.p_old;
      auto& p_now = dcp_param.p_now;
    
      std::vector<T> p_new(b.size());
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] 
	  * e_old + b[i][3] * e_now + b[i][4] * e_new;
      }
      std::copy(p_now.begin(), p_now.end(), p_old.begin());
      std::copy(p_new.begin(), p_new.end(), p_now.begin());
    }
  
  protected:
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<DcpAdeElectricParam<T> > param_list;

  private:
    static const std::string tag; // "DcpAdeElectric"
  }; // template DcpAdeElectric

  template <typename T>
  const std::string DcpAdeElectric<T>::tag = "DcpAdeElectric";

  template <typename T> 
  class DcpAdeEx: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   DcpAdeElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;
      T& e_old = dcp_param.e_old;

      const T& e_now = ex(i,j,k);
      const T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			      (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) 
	+ c[1] * (dps_sum(static_cast<T>(0), dcp_param) + 
		  cps_sum(static_cast<T>(0), dcp_param))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_old, e_now, e_new, dcp_param);
      update_p(e_old, e_now, e_new, dcp_param);
      
      e_old = e_now;
      ex(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::idx_list;
    using DcpAdeElectric<T>::param_list;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  }; // template DcpAdeEx

  template <typename T> 
  class DcpAdeEy: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   DcpAdeElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;
      T& e_old = dcp_param.e_old;
      
      const T& e_now = ey(i,j,k);
      T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * (dps_sum(static_cast<T>(0), dcp_param) + 
		  cps_sum(static_cast<T>(0), dcp_param))
	+ c[2] * e_old + c[3] * e_now;

      update_q(e_old, e_now, e_new, dcp_param);
      update_p(e_old, e_now, e_new, dcp_param);
      
      e_old = e_now;
      ey(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::idx_list;
    using DcpAdeElectric<T>::param_list;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> 
  class DcpAdeEz: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   DcpAdeElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;
      T& e_old = dcp_param.e_old;

      const T& e_now = ez(i,j,k);
      T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * (dps_sum(static_cast<T>(0), dcp_param) + 
		  cps_sum(static_cast<T>(0), dcp_param))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_old, e_now, e_new, dcp_param);
      update_p(e_old, e_now, e_new, dcp_param);
      
      e_old = e_now;
      ez(i,j,k) = e_new;
    }

  protected:
    using DcpAdeElectric<T>::idx_list;
    using DcpAdeElectric<T>::param_list;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  }; // template DcpAdeEz

  template <typename T> 
  class DcpAdeHx: public DielectricHx<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpAdeHx<T>::tag;
    }

  private:
    static const std::string tag; // "DcpMagnetic"
  }; // template DcpAdeHx

  template <typename T>
  const std::string DcpAdeHx<T>::tag = "DcpMagnetic";

  template <typename T> 
  class DcpAdeHy: public DielectricHy<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpAdeHy<T>::tag;
    }

  private:
    static const std::string tag; // "DcpMagnetic"
  }; // template DcpAdeHy

  template <typename T>
  const std::string DcpAdeHy<T>::tag = "DcpMagnetic";

  template <typename T> 
  class DcpAdeHz: public DielectricHz<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpAdeHz<T>::tag;
    }

  private:
    static const std::string tag; // "DcpMagnetic"
  }; // template DcpAdeHz

  template <typename T>
  const std::string DcpAdeHz<T>::tag = "DcpMagnetic";

  /* The following piecewise-linear recursive-convolution(PLRC) 
   * implementation  of the DCP model is based on the following
   * references.
   * 
   * Taflove and S. C. Hagness, Computational Electrodynamics: The
   * Finite-Difference Time-Domain Method, 3rd ed. (Artech House 
   * Publishers, 2005).
   */

  template <typename T> 
  struct DcpPlrcElectricParam: ElectricParam<T>
  {
    std::vector<std::array<double, 3> > a;
    std::vector<std::array<std::complex<double>, 3> > b;
    std::array<double, 3> c;
    // *_re and *_im are for the real and imaginary part of the 
    // e-field, respectively.
    std::vector<double> psi_dp_re, psi_dp_im;
    std::vector<std::complex<double> > psi_cp_re, psi_cp_im;
  }; // template DcpPlrcElectricParam

  template <typename T> 
  struct DcpPlrcMagneticParam: MagneticParam<T>
  {
  }; // template DcpPlrcMagneticParam

  template <typename T> 
  class DcpPlrcElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpPlrcElectric<T>::tag;
    }

    double
    get_eps_inf(const int* const idx, int idx_size) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      const int i = position(index);
      if (i < 0)
	return 0;
      else
	return param_list[i].eps_inf;
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto dcp_param = *static_cast<const DcpPlrcElectricParam<T> * const>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(dcp_param);
      
      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dcp_ptr = static_cast<const DcpPlrcElectric<T>*>(pm_ptr);
      std::copy(dcp_ptr->idx_list.begin(), dcp_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(dcp_ptr->param_list.begin(), dcp_ptr->param_list.end(), std::back_inserter(param_list));
      return this;
    }

    void 
    update_psi_dp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T>& dcp_param)
    {
      const auto& a = dcp_param.a;
      auto& psi_dp_re = dcp_param.psi_dp_re;
      auto& psi_dp_im = dcp_param.psi_dp_im;
      
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	psi_dp_re[i] = a[i][0] * e_new.real() 
	  + a[i][1] * e_now.real() + a[i][2] * psi_dp_re[i];
	psi_dp_im[i] = a[i][0] * e_new.imag() 
	  + a[i][1] * e_now.imag() + a[i][2] * psi_dp_im[i];
      }
    }

    void 
    update_psi_cp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T>& dcp_param)
    {
      const auto& b = dcp_param.b;
      auto& psi_cp_re = dcp_param.psi_cp_re;
      auto& psi_cp_im = dcp_param.psi_cp_im;
      
      for (typename std::vector<std::complex<double> >::size_type i = 0; i < b.size(); ++i) {
	psi_cp_re[i] = b[i][0] * e_new.real() 
	  + b[i][1] * e_now.real() + b[i][2] * psi_cp_re[i];
	psi_cp_im[i] = b[i][0] * e_new.imag()
	  + b[i][1] * e_now.imag() + b[i][2] * psi_cp_im[i];
      }
    }

    std::complex<double> 
    psi_total(const DcpPlrcElectricParam<T>& dcp_param) const
    {
      const auto& psi_dp_re = dcp_param.psi_dp_re;
      const auto& psi_dp_im = dcp_param.psi_dp_im;
      const auto& psi_cp_re = dcp_param.psi_cp_re;
      const auto& psi_cp_im = dcp_param.psi_cp_im;

      double psi_re = 
	std::accumulate(psi_dp_re.begin(), psi_dp_re.end(), 0.0) + 
	std::accumulate(psi_cp_re.begin(), psi_cp_re.end(), std::complex<double>(0)).real();
      
      double psi_im = 
	std::accumulate(psi_dp_im.begin(), psi_dp_im.end(), 0.0) +
	std::accumulate(psi_cp_im.begin(), psi_cp_im.end(), std::complex<double>(0)).real();
      
      return std::complex<double>(psi_re, psi_im);
    }
    
  protected:
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<DcpPlrcElectricParam<T> > param_list;

  private:
    static const std::string tag; // "DcpPlrcElectric"
  }; // template DcpPlrcElectric

  template <typename T>
  const std::string DcpPlrcElectric<T>::tag = "DcpPlrcElectric";

  template <typename S, typename T>
  static inline T& 
  assign(const std::complex<S>& in, T& out)
  {
    return out = static_cast<T>(in.real());
  }

  template <typename S, typename T>
  static inline std::complex<T>& 
  assign(const std::complex<S>& in, std::complex<T>& out)
  {
    return out = in;
  }

  template <typename T> class DcpPlrcEx: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx, 
	   DcpPlrcElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;

      const std::complex<double> e_now = ex(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
		(hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) +
	c[1] * e_now + c[2] * psi_total(dcp_param);
      
      update_psi_dp(e_now, e_new, dcp_param);
      update_psi_cp(e_now, e_new, dcp_param);

      assign(e_new, ex(i,j,k));
  }

  protected:
    using DcpPlrcElectric<T>::idx_list;
    using DcpPlrcElectric<T>::param_list;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  }; // template DcpPlrcEx

  template <typename T> class DcpPlrcEy: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx, 
	   DcpPlrcElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;

      const std::complex<double> e_now = ey(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
		(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx) +
	c[1] * e_now + c[2] * psi_total(dcp_param);

      update_psi_dp(e_now, e_new, dcp_param);
      update_psi_cp(e_now, e_new, dcp_param);

      assign(e_new, ey(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::idx_list;
    using DcpPlrcElectric<T>::param_list;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  }; // template DcpPlrcEy

  template <typename T> class DcpPlrcEz: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n, *idx, *param);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx, 
	   DcpPlrcElectricParam<T>& dcp_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const auto& c = dcp_param.c;

      const std::complex<double> e_now = ez(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
		(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy) +
	c[1] * e_now + c[2] * psi_total(dcp_param);

      update_psi_dp(e_now, e_new, dcp_param);
      update_psi_cp(e_now, e_new, dcp_param);
      
      assign(e_new, ez(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::idx_list;
    using DcpPlrcElectric<T>::param_list;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> 
  class DcpPlrcHx: public DielectricHx<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpPlrcHx<T>::tag;
    }

  private:
    static const std::string tag; // "DcpPlrcMagnetic"
  }; // template DcpPlrcHx

  template <typename T>
  const std::string DcpPlrcHx<T>::tag = "DcpPlrcMagnetic";

  template <typename T> 
  class DcpPlrcHy: public DielectricHy<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpPlrcHy<T>::tag;
    }

  private:
    static const std::string tag; // "DcpPlrcMagnetic"
  }; // template DcpPlrcHy

  template <typename T>
  const std::string DcpPlrcHy<T>::tag = "DcpPlrcMagnetic";

  template <typename T> 
  class DcpPlrcHz: public DielectricHz<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DcpPlrcHz<T>::tag;
    }

  private:
    static const std::string tag; // "DcpPlrcMagnetic"
  }; // template DcpPlrcHz

  template <typename T>
  const std::string DcpPlrcHz<T>::tag = "DcpPlrcMagnetic";
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DCP_HH_*/
