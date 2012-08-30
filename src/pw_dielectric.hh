/* This implementation is based on the following article.
 *
 * K. S. Yee, "Numerical solution of initial boundary value 
 * problems involving Maxwell¡¯s equations in isotropic media," 
 * IEEE Trans. Antennas Propag. 14, 302-307 (1966).
 */

#ifndef PW_DIELECTRIC_HH_
#define PW_DIELECTRIC_HH_

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
  template <typename T> 
  struct DielectricElectricParam: ElectricParam<T>
  {
  }; // template DielectricElectricParam
    
  template <typename T> 
  struct DielectricMagneticParam: MagneticParam<T>
  {
  }; // template DielectricMagneticParam

  template <typename T>
  class DielectricElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DielectricElectric<T>::tag;
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

      const auto& dielectric_param = *static_cast<const DielectricElectricParam<T>*>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(dielectric_param);

      return this;
    }
    
    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dielectric_ptr 
	= static_cast<const DielectricElectric<T>*>(pm_ptr);
      std::copy(dielectric_ptr->idx_list.begin(), 
		dielectric_ptr->idx_list.end(), 
		std::back_inserter(idx_list));
      std::copy(dielectric_ptr->param_list.begin(), 
		dielectric_ptr->param_list.end(), 
		std::back_inserter(param_list));
      return this;
    }

  protected:
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<DielectricElectricParam<T> > param_list;

  private:
    static const std::string tag; // "DielectricElectric"
  }; // template DielectricElectric

  template <typename T>
  const std::string DielectricElectric<T>::tag = "DielectricElectric";

  template <typename T> 
  class DielectricEx: public DielectricElectric<T>
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
	   const DielectricElectricParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double eps_inf = dielectric_param.eps_inf;

      ex(i,j,k) += dt / eps_inf * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
				   (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
    }
    
  protected:
    using DielectricElectric<T>::idx_list;
    using DielectricElectric<T>::param_list;
  }; // template DielectricEx

  template <typename T>
  class DielectricEy: public DielectricElectric<T>
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
	   const DielectricElectricParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double eps_inf = dielectric_param.eps_inf;

      ey(i,j,k) += dt / eps_inf * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
				   (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
    }

  protected:
    using DielectricElectric<T>::idx_list;
    using DielectricElectric<T>::param_list;
  }; // template DielectricEy

  template <typename T> 
  class DielectricEz: public DielectricElectric<T>
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
	       dx, dy, dt, n,
      	       *idx, *param);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n, 
	   const Index3& idx, 
	   const DielectricElectricParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double eps_inf = dielectric_param.eps_inf;

      ez(i,j,k) += dt / eps_inf * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx -
      				   (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
    }

  protected:
    using DielectricElectric<T>::idx_list;
    using DielectricElectric<T>::param_list;
  }; // template DielectricEz

  template <typename T> 
  class DielectricMagnetic: public MaterialMagnetic<T>
  {
  public:
    const std::string& 
    name() const 
    {
      return DielectricMagnetic<T>::tag;
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
      
      const auto& dielectric_param = *static_cast<const DielectricMagneticParam<T>*>(pm_param_ptr);
      
      idx_list.push_back(index);
      param_list.push_back(dielectric_param);
      
      return this;
    }
    
    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dielectric_ptr = static_cast<const DielectricMagnetic<T>*>(pm_ptr);
      std::copy(dielectric_ptr->idx_list.begin(), dielectric_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(dielectric_ptr->param_list.begin(), dielectric_ptr->param_list.end(), std::back_inserter(param_list));
      return this;
    }

  protected:
    using MaterialMagnetic<T>::position;
    using MaterialMagnetic<T>::idx_list;
    std::vector<DielectricMagneticParam<T> > param_list;

  private:
    static const std::string tag; // "DielectricMagnetic"
  }; // template DielectriMagnetic

  template <typename T>
  const std::string DielectricMagnetic<T>::tag = "DielectricMagnetic";

  template <typename T> 
  class DielectricHx: public DielectricMagnetic<T>
  {
  public:
    void
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
	   const DielectricMagneticParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double mu_inf = dielectric_param.mu_inf;

      hx(i,j,k) += dt / mu_inf * ((ey(i,j-1,k) - ey(i,j-1,k-1)) / dz -
      				  (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy);
    }

  protected:
    using DielectricMagnetic<T>::idx_list;
    using DielectricMagnetic<T>::param_list;
  }; // template DielectricHx

  template <typename T> 
  class DielectricHy: public DielectricMagnetic<T>
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

  protected:
    using DielectricMagnetic<T>::idx_list;
    using DielectricMagnetic<T>::param_list;
  }; // template DielectricHy

  template <typename T> 
  class DielectricHz: public DielectricMagnetic<T>
  {
  public:
    void
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
	   const DielectricMagneticParam<T>& dielectric_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      const double mu_inf = dielectric_param.mu_inf;
      
      hz(i,j,k) += dt / mu_inf * ((ex(i-1,j,k) - ex(i-1,j-1,k)) / dy -
				  (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx);
    }

  protected:
    using DielectricMagnetic<T>::idx_list;
    using DielectricMagnetic<T>::param_list;
  }; // template DielectricHz
} // namespace gmes

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_DIELECTRIC_HH_
