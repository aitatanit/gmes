/* This implementation is based on the following article.
 *
 * M. Okoniewski and E. Okoniewska, "Drude dispersion in ADE FDTD revisited,"
 * Electronics Letters, vol. 42, no. 9, pp. 503-504, 2006.
 */

#ifndef PW_DRUDE_HH_
#define PW_DRUDE_HH_

#include <vector>

#include "pw_dielectric.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> struct DrudeElectricParam: public ElectricParam<T>
  {
    std::vector<std::array<double, 3> > a;
    std::array<double, 3> c;
    std::vector<T> q_now, q_new;
  }; // template DrudeElectricParam

  template <typename T> struct DrudeMagneticParam: public MagneticParam<T>
  {
  }; // template DrudeMagneticParam

  template <typename T> class DrudeElectric: public MaterialElectric<T>
  {
  public:
    ~DrudeElectric()
    {
      for (auto v: param) {
	delete static_cast<DrudeElectricParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto drude_param_ptr = static_cast<const DrudeElectricParam<T>* const>(pm_param_ptr);
      auto new_param_ptr = new DrudeElectricParam<T>();
      new_param_ptr->eps_inf = drude_param_ptr->eps_inf;
      std::copy(drude_param_ptr->a.begin(),
		drude_param_ptr->a.end(),
		std::back_inserter(new_param_ptr->a));
      std::copy(drude_param_ptr->c.begin(),
		drude_param_ptr->c.end(),
		new_param_ptr->c.begin());
      new_param_ptr->q_now.resize(drude_param_ptr->a.size(), static_cast<T>(0));
      new_param_ptr->q_new.resize(drude_param_ptr->a.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    };

    T 
    dps_sum(const T& init, const DrudeElectricParam<T>& drude_param) const
    {
      const auto& a = drude_param.a;
      const auto& q_now = drude_param.q_now;
      const auto& q_new = drude_param.q_new;

      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	sum += q_new[i] - q_now[i];
      }

      return sum;
    }

    void 
    update_q(const T& e_now, DrudeElectricParam<T>& drude_param)
    {
      const std::vector<std::array<double, 3> >& a = drude_param.a;
      std::vector<T>& q_now = drude_param.q_now;
      std::vector<T>& q_new = drude_param.q_new;

      std::vector<T> q_old(q_now);
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * e_now;
      }
    }

  protected:
    using MaterialElectric<T>::param;
  }; // template DrudeElectric

  template <typename T> class DrudeEx: public DrudeElectric<T>
  {
  public:
    void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) {
	auto drude_param_ptr = static_cast<DrudeElectricParam<T>*>(v.second);
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n, 
    	       v.first, *drude_param_ptr);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n, 
	   const Index3& idx, 
	   DrudeElectricParam<T>& drude_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = drude_param.c;
      
      const T& e_now = ex(i,j,k);
      update_q(e_now, drude_param);
      ex(i,j,k) = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			  (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
	+ c[1] * dps_sum(static_cast<T>(0), drude_param) + c[2] * e_now;
    }

  protected:
    using DrudeElectric<T>::param;
    using DrudeElectric<T>::update_q;
    using DrudeElectric<T>::dps_sum;
  }; // template DrudeEx

  template <typename T> class DrudeEy: public DrudeElectric<T>
  {
  public:
    void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto drude_param_ptr = static_cast<DrudeElectricParam<T>*>(v.second);
    	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n,
    	       v.first, *drude_param_ptr);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx, 
	   DrudeElectricParam<T>& drude_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = drude_param.c;
      
      const T& e_now = ey(i,j,k);
      update_q(e_now, drude_param);
      ey(i,j,k) = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			  (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * dps_sum(static_cast<T>(0), drude_param) + c[2] * e_now;
    }

  protected:
    using DrudeElectric<T>::param;
    using DrudeElectric<T>::update_q;
    using DrudeElectric<T>::dps_sum;
  }; // template DrudeEy

  template <typename T> class DrudeEz: public DrudeElectric<T>
  {
  public:
    void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto drude_param_ptr = static_cast<DrudeElectricParam<T>*>(v.second);
    	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n,
    	       v.first, *drude_param_ptr);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx, 
	   DrudeElectricParam<T>& drude_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = drude_param.c;
      
      const T& e_now = ez(i,j,k);
      update_q(e_now, drude_param);
      ez(i,j,k) = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			  (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * dps_sum(static_cast<T>(0), drude_param) + c[2] * e_now;
    }

  protected:
    using DrudeElectric<T>::param;
    using DrudeElectric<T>::update_q;
    using DrudeElectric<T>::dps_sum;
  }; // template DrudeEz

  template <typename T> class DrudeHx: public DielectricHx<T>
  {
  }; // template DrudeHx

  template <typename T> class DrudeHy: public DielectricHy<T>
  {
  }; // template DrudeHy

  template <typename T> class DrudeHz: public DielectricHz<T>
  {
  }; // template DrudeHz
} // namespace gmes

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_DRUDE_HH_
