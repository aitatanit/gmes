/* The implementation of Drude-critical points model based on the following
 * articles.
 *
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the 
 * optical properties of gold," The Journal of Chemical Physics, vol. 125,
 * no. 16, pp. 164705-3, Oct. 2006.
 *
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "Erratum: An analytic model
 * for the optical properties of gold" [J. Chem. Phys. 125, 164705 (2006)]," 
 * The Journal of Chemical Physics, vol. 127, no. 18, pp. 189901-1, Nov. 2007.
 *
 * A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-
 * Difference Time-Domain Method, Third Edition, 3rd ed. Artech House 
 * Publishers, 2005.
 */

#ifndef PW_DCP_HH_
#define PW_DCP_HH_

#include <complex>
#include <numeric>
#include <vector>
#include "pw_dielectric.hh"

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  /**************************************************/
  /* Auxiliary Differential Equation Implementation */
  /**************************************************/

  template <typename T> struct DcpAdeElectricParam: ElectricParam<T>
  {
    // parameters for the ADE of the Drude model
    std::vector<std::array<double, 3> > a;
    // parameters for the ADE of critical points model
    std::vector<std::array<double, 4> > b;
    // parameters for the electric field update equations
    std::array<double, 4> c;
    T e_old;
    std::vector<T> q_old, q_now, p_old, p_now;
  };
  
  template <typename T> struct DcpAdeMagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DcpAdeElectric: public MaterialElectric<T>
  {
  public:
    ~DcpAdeElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DcpAdeElectricParam<T> *>(iter->second);
      }
      param.clear();
    }
    
    void 
    attach(const int idx[3], int idx_size,
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<DcpAdeElectricParam<T> *>(iter->second);
	param.erase(iter);
      }

      const DcpAdeElectricParam<T> * const DcpAdeElectricParameter_ptr
	= static_cast<const DcpAdeElectricParam<T> * const>(parameter);
      DcpAdeElectricParam<T> *param_ptr;
      param_ptr = new DcpAdeElectricParam<T>();
      param_ptr->eps = DcpAdeElectricParameter_ptr->eps;
      std::copy(DcpAdeElectricParameter_ptr->a.begin(),
		DcpAdeElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(DcpAdeElectricParameter_ptr->b.begin(),
		DcpAdeElectricParameter_ptr->b.end(),
		std::back_inserter(param_ptr->b));
      std::copy(DcpAdeElectricParameter_ptr->c.begin(),
		DcpAdeElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->e_old = static_cast<T>(0);
      param_ptr->q_old.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->q_now.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->p_old.resize(param_ptr->b.size(), static_cast<T>(0));
      param_ptr->p_now.resize(param_ptr->b.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, param_ptr));
    }
    
    T 
    dps_sum(const T& init, const DcpAdeElectricParam<T> * const ptr) const
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      const std::vector<T>& q_old = ptr->q_old;
      const std::vector<T>& q_now = ptr->q_now;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
      }

      return sum;
    }
      
    T 
    cps_sum(const T& init, const DcpAdeElectricParam<T> * const ptr) const
    {
      const std::vector<std::array<double, 4> >& b = ptr->b;
      const std::vector<T>& p_old = ptr->p_old;
      const std::vector<T>& p_now = ptr->p_now;

      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
      }
      
      return sum;
    }

    void 
    update_q(const T& e_now, DcpAdeElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      std::vector<T>& q_old = ptr->q_old;
      std::vector<T>& q_now = ptr->q_now;

      std::vector<T> q_new(a.size());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * e_now;
      }
      
      std::copy(q_now.begin(), q_now.end(), q_old.begin());
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
    }
    
    void 
    update_p(const T& e_old, const T& e_now, const T& e_new,
	     DcpAdeElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 4> >& b = ptr->b;
      std::vector<T>& p_old = ptr->p_old;
      std::vector<T>& p_now = ptr->p_now;
    
      std::vector<T> p_new(b.size());
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] 
	  * (e_old - e_new) + b[i][3] * e_now;
      }
      std::copy(p_now.begin(), p_now.end(), p_old.begin());
      std::copy(p_new.begin(), p_new.end(), p_now.begin());
    }
  
  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class DcpAdeEx: public DcpAdeElectric<T>
  {
  public:
    void 
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;
 
      const T& e_now = ex(i,j,k);
      const T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			      (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) 
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ex(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeEy: public DcpAdeElectric<T>
  {
  public:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;
      
      const T& e_now = ey(i,j,k);
      T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;

      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ey(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeEz: public DcpAdeElectric<T>
  {
  public:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;

      const T& e_now = ez(i,j,k);
      T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ez(i,j,k) = e_new;
    }

  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeHx: public DielectricHx<T>
  {
  };

  template <typename T> class DcpAdeHy: public DielectricHy<T>
  {
  };

  template <typename T> class DcpAdeHz: public DielectricHz<T>
  {
  };

  /*********************************************************/
  /* Piecewise-Linear Recursive Convolution Implementation */
  /*********************************************************/

  template <typename T> struct DcpPlrcElectricParam: ElectricParam<T>
  {
    std::vector<std::array<double, 3> > a;
    std::vector<std::array<std::complex<double>, 3> > b;
    std::array<double, 3> c;
    // *_re and *_im are for the real and imaginary part of the e-field, 
    // respectively.
    std::vector<double> psi_dp_re, psi_dp_im;
    std::vector<std::complex<double> > psi_cp_re, psi_cp_im;
  };

  template <typename T> struct DcpPlrcMagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DcpPlrcElectric: public MaterialElectric<T>
  {
  public:
    ~DcpPlrcElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DcpPlrcElectricParam<T> *>(iter->second);
      }
      param.clear();
    }

    void
    attach(const int idx[3], int idx_size,
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<DcpPlrcElectricParam<T> *>(iter->second);
	param.erase(iter);
      }

      const DcpPlrcElectricParam<T> * const DcpPlrcElectricParameter_ptr
	= static_cast<const DcpPlrcElectricParam<T> * const>(parameter);
      DcpPlrcElectricParam<T> *param_ptr;
      param_ptr = new DcpPlrcElectricParam<T>();

      param_ptr->eps = DcpPlrcElectricParameter_ptr->eps;
      std::copy(DcpPlrcElectricParameter_ptr->a.begin(),
		DcpPlrcElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(DcpPlrcElectricParameter_ptr->b.begin(),
		DcpPlrcElectricParameter_ptr->b.end(),
		std::back_inserter(param_ptr->b));
      std::copy(DcpPlrcElectricParameter_ptr->c.begin(),
		DcpPlrcElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->psi_dp_re.resize(param_ptr->a.size(), 0);
      param_ptr->psi_dp_im.resize(param_ptr->a.size(), 0);
      param_ptr->psi_cp_re.resize(param_ptr->b.size(), std::complex<double>(0));
      param_ptr->psi_cp_im.resize(param_ptr->b.size(), std::complex<double>(0));
      
      param.insert(std::make_pair(index, param_ptr));
    }

    void 
    update_psi_dp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      std::vector<double>& psi_dp_re = ptr->psi_dp_re;
      std::vector<double>& psi_dp_im = ptr->psi_dp_im;
      
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	psi_dp_re[i] = a[i][0] * e_new.real() + a[i][1] * e_now.real() 
	  + a[i][2] * psi_dp_re[i];
	psi_dp_im[i] = a[i][0] * e_new.imag() + a[i][1] * e_now.imag() 
	  + a[i][2] * psi_dp_im[i];
      }
    }

    void 
    update_psi_cp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T> * const ptr)
    {
      const std::vector<std::array<std::complex<double>, 3> >& b = ptr->b;
      std::vector<std::complex<double> >& psi_cp_re = ptr->psi_cp_re;
      std::vector<std::complex<double> >& psi_cp_im = ptr->psi_cp_im;
      
      for (typename std::vector<std::complex<double> >::size_type i = 0; 
	   i < b.size(); ++i) {
	psi_cp_re[i] = b[i][0] * e_new.real() + b[i][1] * e_now.real()
	  + b[i][2] * psi_cp_re[i];
	psi_cp_im[i] = b[i][0] * e_new.imag() + b[i][1] * e_now.imag()
	  + b[i][2] * psi_cp_im[i];
      }
    }

    std::complex<double> 
    psi_total(const DcpPlrcElectricParam<T> * const ptr) const
    {
      const std::vector<double>& psi_dp_re = ptr->psi_dp_re;
      const std::vector<double>& psi_dp_im = ptr->psi_dp_im;
      const std::vector<std::complex<double> >& psi_cp_re = ptr->psi_cp_re;
      const std::vector<std::complex<double> >& psi_cp_im = ptr->psi_cp_im;

      double psi_re = 
	std::accumulate(psi_dp_re.begin(), psi_dp_re.end(), double(0)) + 
	std::accumulate(psi_cp_re.begin(), psi_cp_re.end(), std::complex<double>(0)).real();
      
      double psi_im = 
	std::accumulate(psi_dp_im.begin(), psi_dp_im.end(), double(0)) +
	std::accumulate(psi_cp_im.begin(), psi_cp_im.end(), std::complex<double>(0)).real();
      
      return std::complex<double>(psi_re, psi_im);
    }
    
  protected:
    using MaterialElectric<T>::param;
  };
  
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
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ex(i,j,k);
      const std::complex<double> e_new = 
	c[0] * e_now + c[1] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			       (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
	+ c[2] * psi_total(ptr);
      
      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);

      assign(e_new, ex(i,j,k));
  }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcEy: public DcpPlrcElectric<T>
  {
  public:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ey(i,j,k);
      const std::complex<double> e_new = 
	c[0] * e_now + c[1] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			       (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[2] * psi_total(ptr);

      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);

      assign(e_new, ey(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcEz: public DcpPlrcElectric<T>
  {
  public:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ez(i,j,k);
      const std::complex<double> e_new = 
	c[0] * e_now + c[1] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			       (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[2] * psi_total(ptr);

      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);
      
      assign(e_new, ez(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcHx: public DielectricHx<T>
  {
  };

  template <typename T> class DcpPlrcHy: public DielectricHy<T>
  {
  };

  template <typename T> class DcpPlrcHz: public DielectricHz<T>
  {
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DCP_HH_*/
