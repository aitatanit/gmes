/* This implementation is based on the following article.
 *
 * S. Gedney, "Perfectly Matched Layer Absorbing Boundary Conditions,"
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * A. Taflove and S.C. Hagness, 3rd ed., Artech House Publishers, 2005,
 * pp. 273-328.
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
  template <typename T> struct UPMLElectricParam: public ElectricParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T d;
  };

  template <typename T> struct UPMLMagneticParam: public MagneticParam<T>
  {
    double c1, c2, c3, c4, c5, c6;
    T b;
  };

  template <typename T> class UPMLElectric: public MaterialElectric<T>
  {
  public:
    ~UPMLElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<UPMLElectricParam<T> *>(iter->second);
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
	delete static_cast<UPMLElectricParam<T> *>(iter->second);
	param.erase(iter);
      }

      const UPMLElectricParam<T>* UPMLElectricParameter_ptr
	= static_cast<const UPMLElectricParam<T> *>(parameter);
      UPMLElectricParam<T> *param_ptr = new UPMLElectricParam<T>();

      param_ptr->eps_inf = UPMLElectricParameter_ptr->eps_inf;
      param_ptr->c1 = UPMLElectricParameter_ptr->c1;
      param_ptr->c2 = UPMLElectricParameter_ptr->c2;
      param_ptr->c3 = UPMLElectricParameter_ptr->c3;
      param_ptr->c4 = UPMLElectricParameter_ptr->c4;
      param_ptr->c5 = UPMLElectricParameter_ptr->c5;
      param_ptr->c6 = UPMLElectricParameter_ptr->c6;
      param_ptr->d = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

  protected:
    using MaterialElectric<T>::param;
  };

  template <typename T> class UPMLEx: public UPMLElectric<T>
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
      
      UPMLElectricParam<T> *ptr
	= static_cast<UPMLElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& d = ptr->d;
      
      const T dstore(d);
      
      d = c1 * d + c2 * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			 (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
      ex(i,j,k) = c3 * ex(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }
    
  protected:
    using UPMLElectric<T>::param;
  };

  template <typename T> class UPMLEy: public UPMLElectric<T>
  {
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const int* const idx, int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UPMLElectricParam<T> *ptr
	= static_cast<UPMLElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& d = ptr->d;
      
      const T dstore(d);

      d = c1 * d + c2 * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			 (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
      ey(i,j,k) = c3 * ey(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }

  protected:
    using UPMLElectric<T>::param;
  };

  template <typename T> class UPMLEz: public UPMLElectric<T>
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
      
      UPMLElectricParam<T> *ptr
	= static_cast<UPMLElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& d = ptr->d;
      
      const T dstore(d);

      d = c1 * d + c2 * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			 (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
      ez(i,j,k) = c3 * ez(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps_inf;
    }

  protected:
    using UPMLElectric<T>::param;
  };

  template <typename T> class UPMLMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~UPMLMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<UPMLMagneticParam<T> *>(iter->second);
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
	delete static_cast<UPMLMagneticParam<T> *>(iter->second);
	param.erase(iter);
      }
      
      const UPMLMagneticParam<T>* UPMLMagneticParameter_ptr 
	= static_cast<const UPMLMagneticParam<T> *>(parameter);
      UPMLMagneticParam<T> *param_ptr = new UPMLMagneticParam<T>();

      param_ptr->mu_inf = UPMLMagneticParameter_ptr->mu_inf;
      param_ptr->c1 = UPMLMagneticParameter_ptr->c1;
      param_ptr->c2 = UPMLMagneticParameter_ptr->c2;
      param_ptr->c3 = UPMLMagneticParameter_ptr->c3;
      param_ptr->c4 = UPMLMagneticParameter_ptr->c4;
      param_ptr->c5 = UPMLMagneticParameter_ptr->c5;
      param_ptr->c6 = UPMLMagneticParameter_ptr->c6;
      param_ptr->b = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class UPMLHx: public UPMLMagnetic<T>
  {
    void 
    update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UPMLMagneticParam<T> *ptr
	= static_cast<UPMLMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& b = ptr->b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy - 
			 (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz);
      hx(i,j,k) = c3 * hx(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UPMLMagnetic<T>::param;
  };

  template <typename T> class UPMLHy: public UPMLMagnetic<T>
  {
  public:
    void 
    update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UPMLMagneticParam<T> *ptr
	= static_cast<UPMLMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& b = ptr->b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz - 
			 (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx);
      hy(i,j,k) = c3 * hy(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UPMLMagnetic<T>::param;
  };

  template <typename T> class UPMLHz: public UPMLMagnetic<T>
  {
  public:
    void 
    update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   double dx, double dy, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UPMLMagneticParam<T> *ptr
	= static_cast<UPMLMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double c1 = ptr->c1;
      double c2 = ptr->c2;
      double c3 = ptr->c3;
      double c4 = ptr->c4;
      double c5 = ptr->c5;
      double c6 = ptr->c6;
      T& b = ptr->b;
      
      const T bstore(b);

      b = c1 * b - c2 * ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx - 
			 (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy);
      hz(i,j,k) = c3 * hz(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu_inf;
    }

  protected:
    using UPMLMagnetic<T>::param;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_UPML_HH_*/
