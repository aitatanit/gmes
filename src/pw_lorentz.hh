/* This implementation is based on the following article.
 *
 * A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-
 * Difference Time-Domain Method, 3rd ed. 685 Canton Street Norwood, MA 02062, 
 * USA: Artech House Publishers, 2005.
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
  };

  template <typename T> struct LorentzMagneticParam: public MagneticParam<T>
  {
  };

  template <typename T> class LorentzElectric: public MaterialElectric<T>
  {
  public:
    ~LorentzElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<LorentzElectricParam<T> *>(iter->second);
      }
      param.clear();
    }
    
    PwMaterial<T> *
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<LorentzElectricParam<T> *>(iter->second);
	param.erase(iter);
      }
      
      const LorentzElectricParam<T> * const LorentzElectricParameter_ptr
	= static_cast<const LorentzElectricParam<T> * const>(parameter);
      LorentzElectricParam<T> *param_ptr;
      param_ptr = new LorentzElectricParam<T>();

      param_ptr->eps_inf = LorentzElectricParameter_ptr->eps_inf;
      std::copy(LorentzElectricParameter_ptr->a.begin(),
		LorentzElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(LorentzElectricParameter_ptr->c.begin(),
		LorentzElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->l_now.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->l_new.resize(param_ptr->a.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, param_ptr));

      return this;
    };

    T 
    lps_sum(const T& init, const LorentzElectricParam<T> * const ptr) const
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      const std::vector<T>& l_now = ptr->l_now;
      const std::vector<T>& l_new = ptr->l_new;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	sum += l_new[i] - l_now[i];
      }

      return sum;
    }

    void 
    update_l(const T& e_now, LorentzElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      std::vector<T>& l_now = ptr->l_now;
      std::vector<T>& l_new = ptr->l_new;

      const std::vector<T> l_old(l_now);
      std::copy(l_new.begin(), l_new.end(), l_now.begin());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i)	{
	l_new[i] = a[i][0] * l_old[i] + a[i][1] * l_now[i] + a[i][2] * e_now;
      }
    }

  protected:
    MaterialElectric<T>::param;
  };

  template <typename T> class LorentzEx: public LorentzElectric<T>
  {
  public:
    void 
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];

      LorentzElectricParam<T> *ptr;
      ptr = static_cast<LorentzElectricParam<T> *>(parameter);
      std::array<double, 3>& c = ptr->c;

      T& e_now = ex(i,j,k);
      update_l(e_now, ptr);
      ex(i,j,k) = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			  (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz)
	+ c[1] * lps_sum(static_cast<T>(0), ptr) + c[2] * e_now;
    }

  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  };

  template <typename T> class LorentzEy: public LorentzElectric<T>
  {
  public:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];

      LorentzElectricParam<T> *ptr;
      ptr = static_cast<LorentzElectricParam<T> *>(parameter);
      std::array<double, 3>& c = ptr->c;
      
      T& e_now = ey(i,j,k);
      update_l(e_now, ptr);
      ey(i,j,k) = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			  (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * lps_sum(static_cast<T>(0), ptr) + c[2] * e_now;
    }
    
  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  };

  template <typename T> class LorentzEz: public LorentzElectric<T>
  {
  public:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];

      LorentzElectricParam<T> *ptr;
      ptr = static_cast<LorentzElectricParam<T> *>(parameter);
      std::array<double, 3>& c = ptr->c;
      
      T& e_now = ez(i,j,k);
      update_l(e_now, ptr);
      ez(i,j,k) = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			  (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * lps_sum(static_cast<T>(0), ptr) + c[2] * e_now;
    }

  protected:
    using LorentzElectric<T>::param;
    using LorentzElectric<T>::update_l;
    using LorentzElectric<T>::lps_sum;
  };

  template <typename T> class LorentzHx: public DielectricHx<T>
  {
  };

  template <typename T> class LorentzHy: public DielectricHy<T>
  {
  };

  template <typename T> class LorentzHz: public DielectricHz<T>
  {
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_LORENTZ_HH_*/
