/* This implementation is based on the following article:
 *
 * M. Okoniewski, M. Mrozowski, and M. A. Stuchly, "Simple 
 * treatment of multi-term dispersion in FDTD," IEEE Microw. 
 * Guided Wave Lett. 7, 121-123 (1997).
 */

#ifndef PW_LORENTZ_HH_
#define PW_LORENTZ_HH_

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
  template <typename T> struct LorentzElectricParam: public ElectricParam<T>
  {
    std::vector<std::array<double, 3> > a;
    std::array<double, 3> c;
    std::vector<T> l_now, l_new;
  }; // template LorentzElectricParam

  template <typename T> struct LorentzMagneticParam: public MagneticParam<T>
  {
  }; // template LorentzMagneticParam

  template <typename T> class LorentzElectric: public MaterialElectric<T>
  {
  public:
    ~LorentzElectric()
    {
      for (auto v: param) {
	delete static_cast<LorentzElectricParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto lorentz_param_ptr = static_cast<const LorentzElectricParam<T> * const>(pm_param_ptr);
      auto new_param_ptr = new LorentzElectricParam<T>();

      new_param_ptr->eps_inf = lorentz_param_ptr->eps_inf;
      std::copy(lorentz_param_ptr->a.begin(),
		lorentz_param_ptr->a.end(),
		std::back_inserter(new_param_ptr->a));
      std::copy(lorentz_param_ptr->c.begin(),
		lorentz_param_ptr->c.end(),
		new_param_ptr->c.begin());
      new_param_ptr->l_now.resize(lorentz_param_ptr->a.size(), static_cast<T>(0));
      new_param_ptr->l_new.resize(lorentz_param_ptr->a.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    };

    T 
    lps_sum(const T& init, const LorentzElectricParam<T>& lorentz_param) const
    {
      const auto& a = lorentz_param.a;
      const std::vector<T>& l_now = lorentz_param.l_now;
      const std::vector<T>& l_new = lorentz_param.l_new;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	sum += l_new[i] - l_now[i];
      }

      return sum;
    }

    void 
    update_l(const T& e_now, LorentzElectricParam<T>& lorentz_param)
    {
      const auto& a = lorentz_param.a;
      std::vector<T>& l_now = lorentz_param.l_now;
      std::vector<T>& l_new = lorentz_param.l_new;

      const std::vector<T> l_old(l_now);
      std::copy(l_new.begin(), l_new.end(), l_now.begin());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	l_new[i] = a[i][0] * l_old[i] + a[i][1] * l_now[i] + a[i][2] * e_now;
      }
    }

  protected:
    MaterialElectric<T>::param;
  }; // template LorentzElectric

  template <typename T> class LorentzEx: public LorentzElectric<T>
  {
  public:
    void
    update_all(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       double dy, double dz, double dt, double n)
    {
      for (auto v: param) {
	auto lorentz_param_ptr = static_cast<LorentzElectricParam<T>*>(v.second);
    	update(ex, ex_x_size, ex_y_size, ex_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       dy, dz, dt, n,
    	       v.first, *lorentz_param_ptr);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx, 
	   LorentzElectricParam<T>& lorentz_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = lorentz_param.c;

      const T& e_now = ex(i,j,k);
      update_l(e_now, lorentz_param);
      ex(i,j,k) = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			  (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
	+ c[1] * lps_sum(static_cast<T>(0), lorentz_param) + c[2] * e_now;
    }

  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  }; // template LorentzEx

  template <typename T> class LorentzEy: public LorentzElectric<T>
  {
  public:
    void
    update_all(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	       double dz, double dx, double dt, double n)
    {
      for (auto v: param) {
	auto lorentz_param_ptr = static_cast<LorentzElectricParam<T>*>(v.second);
	update(ey, ey_x_size, ey_y_size, ey_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       hz, hz_x_size, hz_y_size, hz_z_size,
	       dz, dx, dt, n,
	       v.first, *lorentz_param_ptr);
      }
    }

  private:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx, 
	   LorentzElectricParam<T>& lorentz_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = lorentz_param.c;
      
      const T& e_now = ey(i,j,k);
      update_l(e_now, lorentz_param);
      ey(i,j,k) = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			  (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * lps_sum(static_cast<T>(0), lorentz_param) + c[2] * e_now;
    }
    
  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  }; // template LorentzEy

  template <typename T> class LorentzEz: public LorentzElectric<T>
  {
  public:
    void
    update_all(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	       const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	       const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	       double dx, double dy, double dt, double n)
    {
      for (auto v: param) {
	auto lorentz_param_ptr = static_cast<LorentzElectricParam<T>*>(v.second);
    	update(ez, ez_x_size, ez_y_size, ez_z_size,
	       hy, hy_x_size, hy_y_size, hy_z_size,
	       hx, hx_x_size, hx_y_size, hx_z_size,
	       dx, dy, dt, n,
    	       v.first, *lorentz_param_ptr);
      }
    }

  private:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx, 
	   LorentzElectricParam<T>& lorentz_param)
    {
      const int i = idx[0], j = idx[1], k = idx[2];

      const auto& c = lorentz_param.c;
      
      const T& e_now = ez(i,j,k);
      update_l(e_now, lorentz_param);
      ez(i,j,k) = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			  (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * lps_sum(static_cast<T>(0), lorentz_param) + c[2] * e_now;
    }

  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  }; // template LorentzEz

  template <typename T> class LorentzHx: public DielectricHx<T>
  {
  }; // template LorentzHx

  template <typename T> class LorentzHy: public DielectricHy<T>
  {
  }; // template LorentzHy

  template <typename T> class LorentzHz: public DielectricHz<T>
  {
  }; // template LorentzHz
} // namespace gmes

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_LORENTZ_HH_
