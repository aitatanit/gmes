/* This implementation is based on the following article.
 *
 * K. S. Yee, "Numerical solution of initial boundary value problems involving
 * Maxwell's equations in isotropic media," IEEE Transactions on Antennas and 
 * Propagation, vol. 14, no. 3, pp. 302-307, May. 1966.
 */

#ifndef PW_DIELECTRIC_HH_
#define PW_DIELECTRIC_HH_

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
  template <typename T> struct DielectricElectricParam: ElectricParam<T>
  {
  };
    
  template <typename T> struct DielectricMagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DielectricElectric: public MaterialElectric<T>
  {
  public:
    ~DielectricElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DielectricElectricParam<T>*>(iter->second);
      }
      param.clear();
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<DielectricElectricParam<T>*>(iter->second);
	param.erase(iter);
      }

      DielectricElectricParam<T>* param_ptr;
      param_ptr = new DielectricElectricParam<T>();
      param_ptr->eps_inf 
	= static_cast<const DielectricElectricParam<T>*>(parameter)->eps_inf;

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }
    
  protected:
    using MaterialElectric<T>::param;
  };

  template <typename T> class DielectricEx: public DielectricElectric<T>
  {
  public:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double eps_inf = static_cast<DielectricElectricParam<T>*>(parameter)->eps_inf;

      ex(i,j,k) += dt / eps_inf * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			       (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
    }

  protected:
    using DielectricElectric<T>::param;
  };

  template <typename T> class DielectricEy: public DielectricElectric<T>
  {
  public:
    void 
    update(T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n, 
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double eps_inf = static_cast<DielectricElectricParam<T>*>(parameter)->eps_inf;

      ey(i,j,k) += dt / eps_inf * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			       (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
    }

  protected:
    using DielectricElectric<T>::param;
  };

  template <typename T> class DielectricEz: public DielectricElectric<T>
  {
  public:
    void 
    update(T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n, 
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double eps_inf = static_cast<DielectricElectricParam<T>*>(parameter)->eps_inf;
      
      ez(i,j,k) += dt / eps_inf * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx -
			       (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
    }

  protected:
    using DielectricElectric<T>::param;
  };

  template <typename T> class DielectricMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~DielectricMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DielectricMagneticParam<T>*>(iter->second);
      }
      param.clear();
    }

    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<DielectricMagneticParam<T>*>(iter->second);
	param.erase(iter);
      }

      DielectricMagneticParam<T>* param_ptr;
      param_ptr = new DielectricMagneticParam<T>();
      param_ptr->mu_inf 
	= static_cast<const DielectricMagneticParam<T>*>(parameter)->mu_inf;

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }
    
  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class DielectricHx: public DielectricMagnetic<T>
  {
  public:
    void
    update(T* const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n, 
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double mu_inf = static_cast<DielectricMagneticParam<T>*>(parameter)->mu_inf;

      hx(i,j,k) += dt / mu_inf * ((ey(i,j-1,k) - ey(i,j-1,k-1)) / dz -
			      (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy);
    }

  protected:
    using DielectricMagnetic<T>::param;
  };

  template <typename T> class DielectricHy: public DielectricMagnetic<T>
  {
  public:
    void 
    update(T* const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T* const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n, 
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double mu_inf = static_cast<DielectricMagneticParam<T>*>(parameter)->mu_inf;

      hy(i,j,k) += dt / mu_inf * ((ez(i,j,k-1) - ez(i-1,j,k-1)) / dx -
			      (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz);
    }

  protected:
    using DielectricMagnetic<T>::param;
  };

  template <typename T> class DielectricHz: public DielectricMagnetic<T>
  {
  public:
    void 
    update(T* const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T* const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   double dx, double dy, double dt, double n, 
	   const int* const idx, int idx_size, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double mu_inf = static_cast<DielectricMagneticParam<T>*>(parameter)->mu_inf;
      
      hz(i,j,k) += dt / mu_inf * ((ex(i-1,j,k) - ex(i-1,j-1,k)) / dy -
			      (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx);
    }

  protected:
    using DielectricMagnetic<T>::param;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DIELECTRIC_HH_*/
