/* This implementation is based on the following article:
 *
 * S. D. Gedney, "An anisotropic perfectly matched layer-absorbing
 * medium for the truncation of FDTD lattices," IEEE Trans. 
 * Antennas Propag. 44, 1630-1639 (1996).
 */

#ifndef PW_UPML_HH_
#define PW_UPML_HH_

#include "pw_material.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> 
  struct UpmlElectricParam: public ElectricParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T d;
  }; // template UpmlElectricParam

  template <typename T> 
  struct UpmlMagneticParam: public MagneticParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T b;
  }; // template UpmlMagneticParam

  template <typename T> 
  class UpmlElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return UpmlElectric<T>::tag;
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
      
      const auto& upml_param = *static_cast<const UpmlElectricParam<T>*>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(upml_param);

      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto upml_ptr = static_cast<const UpmlElectric<T>*>(pm_ptr);
      std::copy(upml_ptr->idx_list.begin(), upml_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(upml_ptr->param_list.begin(), upml_ptr->param_list.end(), std::back_inserter(param_list));
      return this;
    }

  protected:
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<UpmlElectricParam<T> > param_list;

  private:
    static const std::string tag; // "UpmlElectric"
  }; // template UpmlElectric

  template <typename T>
  const std::string UpmlElectric<T>::tag = "UpmlElectric";

  template <typename T> 
  class UpmlEx: public UpmlElectric<T>
  {
  public:
    virtual void
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
    using UpmlElectric<T>::idx_list;
    using UpmlElectric<T>::param_list;
  }; // template UpmlEx

  template <typename T> 
  class UpmlEy: public UpmlElectric<T>
  {
    virtual void
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
    using UpmlElectric<T>::idx_list;
    using UpmlElectric<T>::param_list;
  }; // template UpmlEy

  template <typename T> 
  class UpmlEz: public UpmlElectric<T>
  {
  public:
    virtual void
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
    using UpmlElectric<T>::idx_list;
    using UpmlElectric<T>::param_list;
  };

  template <typename T> 
  class UpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    const std::string& 
    name() const
    {
      return UpmlMagnetic<T>::tag;
    }

    double
    get_mu_inf(const int* const idx, int idx_size) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      const int i = position(index);
      if (i < 0)
	return 0;
      else
	return param_list[i].mu_inf;
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      
      const auto& upml_param = *static_cast<const UpmlMagneticParam<T>*>(pm_param_ptr);
      
      idx_list.push_back(index);
      param_list.push_back(upml_param);

      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto upml_ptr = static_cast<const UpmlMagnetic<T>*>(pm_ptr);
      std::copy(upml_ptr->idx_list.begin(), upml_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(upml_ptr->param_list.begin(), upml_ptr->param_list.end(), std::back_inserter(param_list));
      return this;
    }

  protected:
    using MaterialMagnetic<T>::position;
    using MaterialMagnetic<T>::idx_list;
    std::vector<UpmlMagneticParam<T> > param_list;

  private:
    static const std::string tag; // "UpmlMagnetic"
  }; // template UpmlMagnetic

  template <typename T>
  const std::string UpmlMagnetic<T>::tag = "UpmlMagnetic";

  template <typename T> 
  class UpmlHx: public UpmlMagnetic<T>
  {
    virtual void
    update_all(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(hx, hx_x_size, hx_y_size, hx_z_size,
	       ez, ez_x_size, ez_y_size, ez_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       dy, dz, dt, n, *idx, *param);
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
    using UpmlMagnetic<T>::idx_list;
    using UpmlMagnetic<T>::param_list;
  }; // template UpmlHx

  template <typename T> 
  class UpmlHy: public UpmlMagnetic<T>
  {
  public:
    virtual void
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
    using UpmlMagnetic<T>::idx_list;
    using UpmlMagnetic<T>::param_list;
  }; // template UpmlHy

  template <typename T> 
  class UpmlHz: public UpmlMagnetic<T>
  {
  public:
    virtual void
    update_all(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(hz, hz_x_size, hz_y_size, hz_z_size,
	       ey, ey_x_size, ey_y_size, ey_z_size,
	       ex, ex_x_size, ex_y_size, ex_z_size,
	       dx, dy, dt, n, *idx, *param);
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
    using UpmlMagnetic<T>::idx_list;
    using UpmlMagnetic<T>::param_list;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_UPML_HH_*/
