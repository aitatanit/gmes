/* This implementation is based on the following article:
 *
 * S. D. Gedney, "An anisotropic perfectly matched layer-absorbing
 * medium for the truncation of FDTD lattices," IEEE Trans. 
 * Antennas Propag. 44, 1630-1639 (1996).
 */

#ifndef PW_UPML_HH_
#define PW_UPML_HH_

#include <iostream>
#include "pw_material.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> struct UpmlElectricParam: public ElectricParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T d;
  };

  template <typename T> struct UpmlMagneticParam: public MagneticParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T b;
  };

  template <typename T> class UpmlElectric: public MaterialElectric<T>
  {
  public:
    ~UpmlElectric()
    {
      for (auto v: param) {
    	delete static_cast<UpmlElectricParam<T>*>(v.second);
      }
      param.clear();
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto upml_param_ptr = static_cast<const UpmlElectricParam<T>*>(pm_param_ptr);
      auto new_param_ptr = new UpmlElectricParam<T>();

      new_param_ptr->eps_inf = upml_param_ptr->eps_inf;
      new_param_ptr->c1 = upml_param_ptr->c1;
      new_param_ptr->c2 = upml_param_ptr->c2;
      new_param_ptr->c3 = upml_param_ptr->c3;
      new_param_ptr->c4 = upml_param_ptr->c4;
      new_param_ptr->c5 = upml_param_ptr->c5;
      new_param_ptr->c6 = upml_param_ptr->c6;
      new_param_ptr->d = static_cast<T>(0);

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }

  protected:
    using MaterialElectric<T>::param;
  };

  template <typename T> class UpmlEx: public UpmlElectric<T>
  {
  public:
    virtual void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) { 
	auto upml_param_ptr = static_cast<UpmlElectricParam<T>*>(v.second);
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   UpmlElectricParam<T>& upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
     
      const double eps_inf = upml_param.eps_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& d = upml_param.d;
      
      const T dstore(d);
      
      d = c1 * d + c2 * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			 (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
      ex(i,j,k) = c3 * ex(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }
    
  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlEy: public UpmlElectric<T>
  {
    virtual void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto upml_param_ptr = static_cast<UpmlElectricParam<T>*>(v.second);
    	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   UpmlElectricParam<T>& upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const double eps_inf = upml_param.eps_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& d = upml_param.d;
      
      const T dstore(d);

      d = c1 * d + c2 * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			 (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
      ey(i,j,k) = c3 * ey(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }

  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlEz: public UpmlElectric<T>
  {
  public:
    virtual void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto upml_param_ptr = static_cast<UpmlElectricParam<T>*>(v.second);
    	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   UpmlElectricParam<T>& upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const double eps_inf = upml_param.eps_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& d = upml_param.d;
      
      const T dstore(d);

      d = c1 * d + c2 * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			 (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
      ez(i,j,k) = c3 * ez(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }

  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~UpmlMagnetic()
    {
      for (auto v: param) {
    	delete static_cast<UpmlMagneticParam<T>*>(v.second);
      }
      param.clear();
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto upml_param_ptr = static_cast<const UpmlMagneticParam<T>*>(pm_param_ptr);
      auto new_param_ptr = new UpmlMagneticParam<T>();

      new_param_ptr->mu_inf = upml_param_ptr->mu_inf;
      new_param_ptr->c1 = upml_param_ptr->c1;
      new_param_ptr->c2 = upml_param_ptr->c2;
      new_param_ptr->c3 = upml_param_ptr->c3;
      new_param_ptr->c4 = upml_param_ptr->c4;
      new_param_ptr->c5 = upml_param_ptr->c5;
      new_param_ptr->c6 = upml_param_ptr->c6;
      new_param_ptr->b = static_cast<T>(0);

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }

  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class UpmlHx: public UpmlMagnetic<T>
  {
    virtual void
    update_all(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) {
	auto upml_param_ptr = static_cast<UpmlMagneticParam<T>*>(v.second);
    	update(hx, hx_x_size, hx_y_size, hx_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       dy, dz, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx, 
	   UpmlMagneticParam<T> upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const double mu_inf = upml_param.mu_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& b = upml_param.b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy - 
			 (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz);
      hx(i,j,k) = c3 * hx(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UpmlMagnetic<T>::param;
  };

  template <typename T> class UpmlHy: public UpmlMagnetic<T>
  {
  public:
    virtual void
    update_all(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto upml_param_ptr = static_cast<UpmlMagneticParam<T>*>(v.second);
    	update(hy, hy_x_size, hy_y_size, hy_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       dz, dx, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx, 
	   UpmlMagneticParam<T>& upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const double mu_inf = upml_param.mu_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& b = upml_param.b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz - 
			 (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx);
      hy(i,j,k) = c3 * hy(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UpmlMagnetic<T>::param;
  };

  template <typename T> class UpmlHz: public UpmlMagnetic<T>
  {
  public:
    virtual void
    update_all(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto upml_param_ptr = static_cast<UpmlMagneticParam<T>*>(v.second);
    	update(hz, hz_x_size, hz_y_size, hz_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       dx, dy, dt, n,
    	       v.first, *upml_param_ptr);
      }
    }

  private:
    void 
    update(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx, 
	   UpmlMagneticParam<T>& upml_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      
      const double mu_inf = upml_param.mu_inf;
      const double c1 = upml_param.c1;
      const double c2 = upml_param.c2;
      const double c3 = upml_param.c3;
      const double c4 = upml_param.c4;
      const double c5 = upml_param.c5;
      const double c6 = upml_param.c6;
      T& b = upml_param.b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx - 
			 (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy);
      hz(i,j,k) = c3 * hz(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UpmlMagnetic<T>::param;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_UPML_HH_*/
