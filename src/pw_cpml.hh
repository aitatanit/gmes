/* This implementation is based on the following article.
 *
 * S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,"
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * A. Taflove and S.C. Hagness, 3rd eds., Artech House Publishers, 2005, 
 * pp. 273-328.
 */ 

#ifndef PW_CPML_HH_
#define PW_CPML_HH_

#include <iostream>
#include <utility>
#include "pw_material.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> struct CpmlElectricParam: public ElectricParam<T>
  {
    double b1, b2, c1, c2, kappa1, kappa2;
    T psi1, psi2;
  }; // template CpmlElectricParam
  
  template <typename T> struct CpmlMagneticParam: public MagneticParam<T>
  {
    double b1, b2, c1, c2, kappa1, kappa2;
    T psi1, psi2;
  }; // template CpmlMagneticParam
  
  template <typename T> class CpmlElectric: public MaterialElectric<T>
  {
  public:
    ~CpmlElectric()
    {
      for (auto v: param) {
	delete static_cast<CpmlElectricParam<T>*>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto cpml_param_ptr
	= static_cast<const CpmlElectricParam<T>* const>(pm_param_ptr);
      auto new_param_ptr = new CpmlElectricParam<T>();
      new_param_ptr->eps_inf = cpml_param_ptr->eps_inf;
      new_param_ptr->b1 = cpml_param_ptr->b1;
      new_param_ptr->b2 = cpml_param_ptr->b2;
      new_param_ptr->c1 = cpml_param_ptr->c1;
      new_param_ptr->c2 = cpml_param_ptr->c2;
      new_param_ptr->kappa1 = cpml_param_ptr->kappa1;
      new_param_ptr->kappa2 = cpml_param_ptr->kappa2;
      new_param_ptr->psi1 = static_cast<T>(0);
      new_param_ptr->psi2 = static_cast<T>(0);

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }

  protected:
    using PwMaterial<T>::param;
  }; // template CpmlElectric

  template <typename T> class CpmlEx: public CpmlElectric<T>
  {
  public:
    void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlElectricParam<T>*>(v.second);
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   CpmlElectricParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double eps_inf = cpml_param.eps_inf;
      const double by = cpml_param.b1;
      const double bz = cpml_param.b2;
      const double cy = cpml_param.c1;
      const double cz = cpml_param.c2;
      const double kappay = cpml_param.kappa1;
      const double kappaz = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = by * psi1 + cy * (hz(i+1,j+1,k) - hz(i+1,j,k)) / dy;
      psi2 = bz * psi2 + cz * (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;
      
      ex(i,j,k) += dt / eps_inf * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy / kappay -
				   (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz / kappaz +
				   psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  }; // template CpmlEx

  template <typename T> class CpmlEy: public CpmlElectric<T>
  {
  public:
    void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlElectricParam<T>*>(v.second);
    	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   CpmlElectricParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double eps_inf = cpml_param.eps_inf;
      const double bz = cpml_param.b1;
      const double bx = cpml_param.b2;
      const double cz = cpml_param.c1;
      const double cx = cpml_param.c2;
      const double kappaz = cpml_param.kappa1;
      const double kappax = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = bz * psi1 + cz * (hx(i,j+1,k+1) - hx(i,j+1,k)) / dz;
      psi2 = bx * psi2 + cx * (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx;
      
      ey(i,j,k) += dt / eps_inf * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz / kappaz -
				   (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx / kappax +
				   psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  }; // template CpmlEy

  template <typename T> class CpmlEz: public CpmlElectric<T>
  {
  public:
    void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlElectricParam<T>*>(v.second);
	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }
    
  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   CpmlElectricParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double eps_inf = cpml_param.eps_inf;
      const double bx = cpml_param.b1;
      const double by = cpml_param.b2;
      const double cx = cpml_param.c1;
      const double cy = cpml_param.c2;
      const double kappax = cpml_param.kappa1;
      const double kappay = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = bx * psi1 + cx * (hy(i+1,j,k+1) - hy(i,j,k+1)) / dx;
      psi2 = by * psi2 + cy * (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy;
      
      ez(i,j,k) += dt / eps_inf * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx / kappax -
				   (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy / kappay +
				   psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  }; // template CpmlEz

  template <typename T> class CpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~CpmlMagnetic()
    {
      for (auto v: param) {
	delete static_cast<CpmlMagneticParam<T>*>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto cpml_param_ptr = static_cast<const CpmlMagneticParam<T>* const>(pm_param_ptr);
      auto new_param_ptr = new CpmlMagneticParam<T>();
      new_param_ptr->mu_inf = cpml_param_ptr->mu_inf;
      new_param_ptr->b1 = cpml_param_ptr->b1;
      new_param_ptr->b2 = cpml_param_ptr->b2;
      new_param_ptr->c1 = cpml_param_ptr->c1;
      new_param_ptr->c2 = cpml_param_ptr->c2;
      new_param_ptr->kappa1 = cpml_param_ptr->kappa1;
      new_param_ptr->kappa2 = cpml_param_ptr->kappa2;
      new_param_ptr->psi1 = static_cast<T>(0);
      new_param_ptr->psi2 = static_cast<T>(0);

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }

  protected:
    using PwMaterial<T>::param;
  }; // template CpmlMagnetic

  template <typename T> class CpmlHx: public CpmlMagnetic<T>
  {
  public:
    void
    update_all(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlMagneticParam<T>*>(v.second);
    	update(hx, hx_x_size, hx_y_size, hx_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       dy, dz, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }

  private:
    void 
    update(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   CpmlMagneticParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double mu_inf = cpml_param.mu_inf;
      const double by = cpml_param.b1;
      const double bz = cpml_param.b2;
      const double cy = cpml_param.c1;
      const double cz = cpml_param.c2;
      const double kappay = cpml_param.kappa1;
      const double kappaz = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = by * psi1 + cy * (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy;
      psi2 = bz * psi2 + cz * (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz;

      hx(i,j,k) -= dt / mu_inf * ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy / kappay -
				  (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz / kappaz +
				  psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::param;
  }; // template CpmlHx

  template <typename T> class CpmlHy: public CpmlMagnetic<T>
  {
  public:
    void
    update_all(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlMagneticParam<T>*>(v.second);
    	update(hy, hy_x_size, hy_y_size, hy_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       dz, dx, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }

  private:
    void 
    update(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   CpmlMagneticParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double mu_inf = cpml_param.mu_inf;
      const double bz = cpml_param.b1;
      const double bx = cpml_param.b2;
      const double cz = cpml_param.c1;
      const double cx = cpml_param.c2;
      const double kappaz = cpml_param.kappa1;
      const double kappax = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = bz * psi1 + cz * (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz;
      psi2 = bx * psi2 + cx * (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx;

      hy(i,j,k) -= dt / mu_inf * ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz / kappaz -
				  (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx / kappax +
				  psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::param;
  }; // template CpmlHy

  template <typename T> class CpmlHz: public CpmlMagnetic<T>
  {
  public:
    void
    update_all(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto cpml_param_ptr = static_cast<CpmlMagneticParam<T>*>(v.second);
    	update(hz, hz_x_size, hz_y_size, hz_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       dx, dy, dt, n,
    	       v.first, *cpml_param_ptr);
      }
    }

  protected:
    using CpmlMagnetic<T>::param;

  private:
    void 
    update(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   CpmlMagneticParam<T>& cpml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const double mu_inf = cpml_param.mu_inf;
      const double bx = cpml_param.b1;
      const double by = cpml_param.b2;
      const double cx = cpml_param.c1;
      const double cy = cpml_param.c2;
      const double kappax = cpml_param.kappa1;
      const double kappay = cpml_param.kappa2;
      T& psi1 = cpml_param.psi1;
      T& psi2 = cpml_param.psi2;

      psi1 = bx * psi1 + cx * (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx;
      psi2 = by * psi2 + cy * (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy;
      
      hz(i,j,k) -= dt / mu_inf * ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx / kappax -
				  (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy / kappay +
				  psi1 - psi2);
    }
  }; // template CpmlHz
} // namespace gmes

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_CPML_HH_
