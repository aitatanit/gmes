/* This implementation is based on the following articles:
 *
 * 1. Ziolkowski, R. W., Arnold, J. M. & Gogny, D. M. 
 * Ultrafast pulse interactions with two-level atoms. 
 * Phys. Rev. A 52, 3082–3094 (1995).
 *
 * 2. Schlottau, F., Piket-May, M. & Wagner, K. Modeling 
 * of femtosecond pulse interaction with inhomogeneously 
 * broadened media using an iterative predictor corrector 
 * FDTD method. Opt. Express 13, 182–194 (2005).
 *
 * This module just handles 1D case with Ex and Hy fields.
 *
 * TODOs
 * 1. Replace array with valarray.
 * 2. 3-level medium
 */

#ifndef PW_DM2_HH_
#define PW_DM2_HH_

#include <array>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <iostream>

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
  // TODO: l2_norm should be a local function.
  template <typename T>
  double 
  l2_norm(const T& e, const std::vector<std::array<T, 3> >& u)
  {
    T accu = e * e;
    for (std::size_t i = 0; i < u.size(); ++i) {  
      accu += std::inner_product(u[i].begin(), u[i].end(), 
                                 u[i].begin(), static_cast<T>(0));
    }

    // accu += std::accumulate(u.begin(), u.end(), 0,
    //                         [](int a, const std::array<int, 3>& b)
    //                         { 
    //                           return a + std::inner_product(b.begin(), b.end(), 
    //                                                         b.begin(),
    //                                                         static_cast<T>(0));
    //                         }
    //                         );

    return std::sqrt(accu);
  } // template l2_norm


  // TODO: rel_error should be a local function.
  template <typename T>
  double
  rel_error(const T& e_new, const std::vector<std::array<T, 3> >& u_new,
            const T& e_ref, const std::vector<std::array<T, 3> >& u_ref)
  {
    T e_diff = e_new - e_ref;

    std::vector<std::array<T, 3> > u_diff;
    u_diff.reserve(u_ref.size());

    for (std::size_t  i = 0; i < u_ref.size(); ++i) {
      std::set_difference(u_new[i].begin(), u_new[i].end(),
                          u_ref[i].begin(), u_ref[i].end(),
                          u_diff[i].begin());
    }

    double denom = l2_norm(e_diff, u_diff);
    double num = l2_norm(e_ref, u_ref);
    
    return denom / num;
  } // template rel_error


  template <typename T> 
  struct Dm2ElectricParam: public ElectricParam<T>
  {
    std::vector<double> omega;
    std::vector<double> n_atom;
    double rho30;
    double gamma;
    double t1, t2;
    double hbar;
    double rtol;

    std::vector<std::array<T, 3> > u;
  }; // template Dm2ElectricParam
  

  template <typename T>
  struct Dm2MagneticParam: public MagneticParam<T>
  {
  }; // template Dm2MagneticParam

  
  template <typename T> 
  class Dm2Electric: public MaterialElectric<T>
  {
  public:
    // TODO: range check of bin
    T
    get_rho(const int* const idx, int idx_size, int bin, int rho_idx, double t) const
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
            return exp(-t / p.t2) * p.u[bin][rho_idx];
	  case 2:
            return p.rho30 + exp(-t / p.t1) * p.u[bin][rho_idx];
	  default:
            throw std::out_of_range("rho_idx should be in [0, 2]");
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

    void
    a(double t, const Dm2ElectricParam<T>& dm2_param,
      std::vector<double>& a_out) const 
    {
      const auto& n_atom = dm2_param.n_atom;
      const auto& gamma = dm2_param.gamma;
      const auto& t2 = dm2_param.t2;

      for (std::size_t i = 0; i < n_atom.size(); ++i) {
        a_out.push_back(n_atom[i] * gamma / t2 * exp(-t / t2));
      }
    }

    void
    b(double t, const Dm2ElectricParam<T>& dm2_param,
      std::vector<double>& b_out) const
    {
      const auto& n_atom = dm2_param.n_atom;
      const auto& gamma = dm2_param.gamma;
      const auto& omega = dm2_param.omega;
      const auto& t2 = dm2_param.t2;

      for (std::size_t i = 0; i < n_atom.size(); ++i) {
        b_out.push_back(n_atom[i] * gamma * omega[i] * exp(-t / t2));
      }
    }

    void 
    c_plus(double t, const Dm2ElectricParam<T>& dm2_param,
           double& c_out) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& t1 = dm2_param.t1;
      const auto& t2 = dm2_param.t2;
      const auto& hbar = dm2_param.hbar;
      
      c_out = 2 / hbar * gamma * exp(-t * (1 / t1 - 1 / t2));
    }

    void
    c_minus(double t, const Dm2ElectricParam<T>& dm2_param,
            double& c_out) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& t1 = dm2_param.t1;
      const auto& t2 = dm2_param.t2;
      const auto& hbar = dm2_param.hbar;

      c_out = 2 / hbar * gamma * exp(-t * (1 / t2 - 1 / t1));
    }

    void
    d(double t, const Dm2ElectricParam<T>& dm2_param,
      double& d_out) const
    {
      const auto& gamma = dm2_param.gamma;
      const auto& rho30 = dm2_param.rho30;
      const auto& t2 = dm2_param.t2;
      const auto& hbar = dm2_param.hbar;

      d_out = 2 / hbar * gamma * rho30 * exp(t / t2);
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
      const std::vector<double>& omega = dm2_param.omega;
      const double rtol = dm2_param.rtol;
      std::vector<std::array<T, 3> >& u = dm2_param.u;
      
      const double t = (n + 0.5) * dt;
      
      std::vector<double> a, b;
      double c_plus, c_minus, d;

      this->a(t, dm2_param, a);
      this->b(t, dm2_param, b);
      this->c_plus(t, dm2_param, c_plus);
      this->c_minus(t, dm2_param, c_minus);
      this->d(t, dm2_param, d);
      
      T e_new = ex(i,j,k);
      std::vector<std::array<T, 3> > u_new = u;

      const T e_old = ex(i,j,k);
      const T hy_dz = (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;

      double error;
      do {
        T e_tmp = e_new;
        std::vector<std::array<T, 3> > u_tmp = u_new;
        
	e_new = e_old - dt * hy_dz;
        for (std::size_t i = 0; i < u.size(); ++i) {
	  e_new -= .5 * dt * a[i] * (u_new[i][0] + u[i][0]);
	  e_new += .5 * dt * b[i] * (u_new[i][1] + u[i][1]);
        }

        for (std::size_t i = 0; i < u.size(); ++i) {
          u_new[i][0] = u[i][0] 
            + .5 * dt * omega[i] * (u_new[i][1] + u[i][1]);
          u_new[i][1] = u[i][1] 
            - .5 * dt * omega[i] * (u_new[i][0] + u[i][0])
            + .25 * dt * c_plus * (u_new[i][2] + u[i][2]) * (e_new + e_old) 
            + .5 * dt * d * (e_new + e_old);
          u_new[i][2] = u[i][2] 
            - .25 * dt * c_minus * (u_new[i][1] + u[i][1]) * (e_new + e_old);
        }
        
        error = rel_error(e_new, u_new, e_tmp, u_tmp);
      } while (error > rtol);
      
      ex(i,j,k) = e_new;
      std::copy(u_new.begin(), u_new.end(), u.begin());
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
