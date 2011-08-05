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

#define ex(i,j,k) ex[((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[((i)*hz_y_size+(j))*hz_z_size+(k)]

namespace gmes
{
  template <typename T> struct UpmlElectricParam: public ElectricParam 
  {
    double c1, c2, c3, c4, c5, c6;
    T d;
  };

  template <typename T> struct UpmlMagneticParam: public MagneticParam 
  {
    double c1, c2, c3, c4, c5, c6;
    T b;
  };

  template <typename T> class UpmlElectric: public MaterialElectric<T>
  {
  public:
    ~UpmlElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<UpmlElectricParam<T> *>(iter->second);
      }
      param.clear();
    }

    void 
    attach(const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<UpmlElectricParam<T> *>(iter->second);
	param.erase(iter);
      }

      const UpmlElectricParam<T>* UpmlElectricParameter_ptr
	= static_cast<const UpmlElectricParam<T> *>(parameter);
      UpmlElectricParam<T> *param_ptr = new UpmlElectricParam<T>();
      param_ptr->eps = UpmlElectricParameter_ptr->eps;
      param_ptr->c1 = UpmlElectricParameter_ptr->c1;
      param_ptr->c2 = UpmlElectricParameter_ptr->c2;
      param_ptr->c3 = UpmlElectricParameter_ptr->c3;
      param_ptr->c4 = UpmlElectricParameter_ptr->c4;
      param_ptr->c5 = UpmlElectricParameter_ptr->c5;
      param_ptr->c6 = UpmlElectricParameter_ptr->c6;
      param_ptr->d = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));
    }

  protected:
    using MaterialElectric<T>::param;
  };

  template <typename T> class UpmlEx: public UpmlElectric<T>
  {
  public:
    void 
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlElectricParam<T> *ptr
	= static_cast<UpmlElectricParam<T> *>(parameter);
      double eps = ptr->eps;
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
      ex(i,j,k) = c3 * ex(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps;
    }
    
  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlEy: public UpmlElectric<T>
  {
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlElectricParam<T> *ptr
	= static_cast<UpmlElectricParam<T> *>(parameter);
      double eps = ptr->eps;
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
      ey(i,j,k) = c3 * ey(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps;
    }

  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlEz: public UpmlElectric<T>
  {
  public:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlElectricParam<T> *ptr
	= static_cast<UpmlElectricParam<T> *>(parameter);
      double eps = ptr->eps;
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
      ez(i,j,k) = c3 * ez(i,j,k) + c4 * (c5 * d - c6 * dstore) / eps;
    }

  protected:
    using UpmlElectric<T>::param;
  };

  template <typename T> class UpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~UpmlMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<UpmlMagneticParam<T> *>(iter->second);
      }
      param.clear();
    }

    void
    attach(const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      MapType::const_iterator iter = param.find(index);
      if (iter != param.end()) {
	std::cerr << "Overwriting the existing index." << std::endl;
	delete static_cast<UpmlMagneticParam<T> *>(iter->second);
	param.erase(iter);
      }
      
      const UpmlMagneticParam<T>* UpmlMagneticParameter_ptr 
	= static_cast<const UpmlMagneticParam<T> *>(parameter);
      UpmlMagneticParam<T> *param_ptr = new UpmlMagneticParam<T>();
      param_ptr->mu = UpmlMagneticParameter_ptr->mu;
      param_ptr->c1 = UpmlMagneticParameter_ptr->c1;
      param_ptr->c2 = UpmlMagneticParameter_ptr->c2;
      param_ptr->c3 = UpmlMagneticParameter_ptr->c3;
      param_ptr->c4 = UpmlMagneticParameter_ptr->c4;
      param_ptr->c5 = UpmlMagneticParameter_ptr->c5;
      param_ptr->c6 = UpmlMagneticParameter_ptr->c6;
      param_ptr->b = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));
    }

  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class UpmlHx: public UpmlMagnetic<T>
  {
    void 
    update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlMagneticParam<T> *ptr
	= static_cast<UpmlMagneticParam<T> *>(parameter);
      double mu = ptr->mu;
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
      hx(i,j,k) = c3 * hx(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
    }

  protected:
    using UpmlMagnetic<T>::param;
  };

  template <typename T> class UpmlHy: public UpmlMagnetic<T>
  {
  public:
    void 
    update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   double dz, double dx, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlMagneticParam<T> *ptr
	= static_cast<UpmlMagneticParam<T> *>(parameter);
      double mu = ptr->mu;
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
      hy(i,j,k) = c3 * hy(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
    }

  protected:
    using UpmlMagnetic<T>::param;
  };

  template <typename T> class UpmlHz: public UpmlMagnetic<T>
  {
  public:
    void 
    update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   double dx, double dy, double dt, double n,
	   const int idx[3], int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      UpmlMagneticParam<T> *ptr
	= static_cast<UpmlMagneticParam<T> *>(parameter);
      double mu = ptr->mu;
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
      hz(i,j,k) = c3 * hz(i,j,k) + c4 * (c5 * b - c6 * bstore) / mu;
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
