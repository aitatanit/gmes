/* This implementation is based on the following article:
 *
 * Ziolkowski, R. W., Arnold, J. M. & Gogny, D. M. 
 * "Ultrafast pulse interactions with two-level atoms. 
 * Phys. Rev. A 52, 3082-3094 (1995).
 *
 * This module just handles 1D case with Ex and Hy fields.
 */

#ifndef PW_DM2_HH_
#define PW_DM2_HH_

#include <array>
#include <vector>
#include <cmath>
#include "pw_dielectric.hh"

using namespace std;

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> 
  struct Dm2ElectricParam: public ElectricParam<T>
  {
    double omega0;
    double rho30;
    double n_atom;
    double gamma;
    double t1, t2;

    std::array<T, 3> u;
  }; // template Dm2ElectricParam
  

  template <typename T> 
  struct Dm2MagneticParam: public MagneticParam<T>
  {
  }; // template Dm2MagneticParam

  
  template <typename T> 
  class Dm2Electric: public MaterialElectric<T>
  {
  public:
    T
    get_rho(const int* const idx, int idx_size, int rho_idx, double t) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      const int i = position(index);

      if (i < 0)
	return 0;
      else {
	const auto& p = param_list[i];
	  
	switch (rho_idx)
	  {
	  case 0:
	  case 1:
	    return exp(-t / p.t2) * p.u[rho_idx];
	  case 2:
	    return p.rho30 + exp(-t / p.t1) * p.u[rho_idx];
	  default:
	    return 0;
	  }
      }
    }
    
    const std::string& 
    name() const
    {
      return Dm2Electric<T>::tag;
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

      const auto& dm2_param = *static_cast<const Dm2ElectricParam<T> * const>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(dm2_param);

      return this;
    };

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dm2_ptr 
	= static_cast<const Dm2Electric<T>*>(pm_ptr);
      std::copy(dm2_ptr->idx_list.begin(), 
		dm2_ptr->idx_list.end(), 
		std::back_inserter(idx_list));
      std::copy(dm2_ptr->param_list.begin(), 
		dm2_ptr->param_list.end(), 
		std::back_inserter(param_list));
      return this;
    }
    
  protected:
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<Dm2ElectricParam<T> > param_list;

    double
    a(double t, const Dm2ElectricParam<T>& dm2_param) const 
    {
      const auto& n_atom = dm2_param.n_atom;
      const auto& gamma = dm2_param.gamma;
      const auto& t2 = dm2_param.t2;

      return n_atom * gamma / t2 * exp(-t / t2);
    }

    double
    b(double t, const Dm2ElectricParam<T>& dm2_param) const
    {
      const auto& n_atom = dm2_param.n_atom;
      const auto& gamma = dm2_param.gamma;
      const auto& omega0 = dm2_param.omega0;
      const auto& t2 = dm2_param.t2;

      return n_atom * gamma * omega0 * exp(-t / t2);
    }

    double 
    c_plus(double t, const Dm2ElectricParam<T>& dm2_param) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& t1 = dm2_param.t1;
      const auto& t2 = dm2_param.t2;

      return 2 * gamma * exp(-t * (1 / t1 - 1 / t2));
    }

    double
    c_minus(double t, const Dm2ElectricParam<T>& dm2_param) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& t1 = dm2_param.t1;
      const auto& t2 = dm2_param.t2;

      return 2 * gamma * exp(- t * (1 / t2 - 1 / t1));
    }

    double
    d(double t, const Dm2ElectricParam<T>& dm2_param) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& rho30 = dm2_param.rho30;
      const auto& t2 = dm2_param.t2;

      return 2 * gamma * rho30 * exp(t / t2);
    }
    
  private:
    static const std::string tag; // "Dm2Electric"
  }; // template Dm2Electric

  template <typename T>
  const std::string Dm2Electric<T>::tag = "Dm2Electric";

  
  template <typename T> 
  class Dm2Ex: public Dm2Electric<T>
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
	   Dm2ElectricParam<T>& dm2_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double omega0 = dm2_param.omega0;
      std::array<T, 3>& u = dm2_param.u;
      
      const double t = (n + 0.5) * dt;
      const double a = this->a(t, dm2_param);
      const double b = this->b(t, dm2_param);
      const double c_plus = this->c_plus(t, dm2_param);
      const double c_minus = this->c_minus(t, dm2_param);
      const double d = this->d(t, dm2_param);
      
      std::array<T, 4> u_new, u_tmp;
      u_new[0] = ex(i,j,k);
      std::copy(u.begin(), u.end(), next(u_new.begin()));

      const T e_old = ex(i,j,k);
      const T hy_dz = (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;

      do {
        std::copy(u_new.begin(), u_new.end(), u_tmp.begin());
        
	u_new[0] = e_old - dt * hy_dz
	  - .5 * dt * a * (u_new[1] + u[0])
	  + .5 * dt * b * (u_new[2] + u[1]);

	u_new[1] = u[0] 
	  + .5 * dt * omega0 * (u_new[2] + u[1]);
	u_new[2] = u[1] 
	  - .5 * dt * omega0 * (u_new[1] + u[0])
	  + .25 * dt * c_plus * (u_new[3] + u[2]) * (u_new[0] + e_old) 
	  + .5 * dt * d * (u_new[0] + e_old);
	u_new[3] = u[2] 
	  - .25 * dt * c_minus * (u_new[2] + u[1]) * (u_new[0] + e_old);
        
      // Ziolkowski tests the saturation with norm of u. However,
      // I tests each element of u, for the sake of simplicity.
      } while(!std::equal(u_new.begin(), u_new.end(), u_tmp.begin(), 
                          [] (const T& a, const T& b) 
                          { return std::abs(a - b) <= 1e-5 * std::abs(b); }));
      
      ex(i,j,k) = u_new[0];
      std::copy(next(u_new.begin()), u_new.end(), u.begin());
    }
    
  protected:
    using Dm2Electric<T>::idx_list;
    using Dm2Electric<T>::param_list;
  }; // template Dm2Ex
  

  template <typename T>
  class Dm2Ey: public Dm2Electric<T>
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
	   const Dm2ElectricParam<T>& dm2_param) const
    {
      // const int i = idx[0], j = idx[1], k = idx[2];

      // ey(i,j,k) += dt * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
      // 			 (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
    }

  protected:
    using Dm2Electric<T>::idx_list;
    using Dm2Electric<T>::param_list;
  }; // template Dm2Ey


  template <typename T> 
  class Dm2Ez: public Dm2Electric<T>
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
	   const Dm2ElectricParam<T>& dm2_param) const
    {
      // const int i = idx[0], j = idx[1], k = idx[2];

      // ez(i,j,k) += dt * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx -
      // 			 (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
    }

  protected:
    using Dm2Electric<T>::idx_list;
    using Dm2Electric<T>::param_list;
  }; // template Dm2Ez


  template <typename T>
  class Dm2Hx: public DielectricHx<T>
  {
  public:
    const std::string&
    name() const
    {
      return Dm2Hx<T>::tag;
    }

  private:
    static const std::string tag; // "Dm2Magnetic"
  }; // template Dm2Hx

  template <typename T>
  const std::string Dm2Hx<T>::tag = "Dm2Magnetic";


  template <typename T> 
  class Dm2Hy: public DielectricHy<T>
  {
  public:
    void
    update_all(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
      	update(hy, hy_x_size, hy_y_size, hy_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       dz, dx, dt, n, *idx, *param);
      }
    }

    const std::string& 
    name() const
    {
      return Dm2Hy<T>::tag;
    }

  private:
    void 
    update(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n, 
	   const Index3& idx, 
	   const DielectricMagneticParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double mu_inf = dielectric_param.mu_inf;
      
      hy(i,j,k) += dt / mu_inf * ((ez(i,j,k-1) - ez(i-1,j,k-1)) / dx -
      				  (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz);
    }
    
    static const std::string tag; // "Dm2Magnetic"

  protected:
    using DielectricHy<T>::idx_list;
    using DielectricHy<T>::param_list;
  }; // template Dm2Hy

  template <typename T>
  const std::string Dm2Hy<T>::tag = "Dm2Magnetic";


  template <typename T>
  class Dm2Hz: public DielectricHz<T>
  {
  public:
    const std::string& 
    name() const
    {
      return Dm2Hz<T>::tag;
    }

  private:
    static const std::string tag; // "Dm2Magnetic"
  }; // template Dm2Hz

  template <typename T>
  const std::string Dm2Hz<T>::tag = "Dm2Magnetic";
} // namespace gmes

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_DM2_HH_
