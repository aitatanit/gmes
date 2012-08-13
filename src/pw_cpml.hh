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
  };
  
  template <typename T> struct CpmlMagneticParam: public MagneticParam<T>
  {
    double b1, b2, c1, c2, kappa1, kappa2;
    T psi1, psi2;
  };
  
  template <typename T> class CpmlElectric: public MaterialElectric<T>
  {
  public:
    ~CpmlElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<CpmlElectricParam<T> *>(iter->second);
      }
      param.clear();
    }
    
    PwMaterial<T> *
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::iterator iter = param.find(index);
      if (iter != param.end()) {
      	std::cerr << "Overwriting the existing index." << std::endl;
      	delete static_cast<CpmlElectricParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const CpmlElectricParam<T> * const CpmlElectricParameter_ptr
	= static_cast<const CpmlElectricParam<T> * const>(parameter);
      CpmlElectricParam<T> *param_ptr;
      param_ptr = new CpmlElectricParam<T>();
      param_ptr->eps_inf = CpmlElectricParameter_ptr->eps_inf;
      param_ptr->b1 = CpmlElectricParameter_ptr->b1;
      param_ptr->b2 = CpmlElectricParameter_ptr->b2;
      param_ptr->c1 = CpmlElectricParameter_ptr->c1;
      param_ptr->c2 = CpmlElectricParameter_ptr->c2;
      param_ptr->kappa1 = CpmlElectricParameter_ptr->kappa1;
      param_ptr->kappa2 = CpmlElectricParameter_ptr->kappa2;
      param_ptr->psi1 = static_cast<T>(0);
      param_ptr->psi2 = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class CpmlEx: public CpmlElectric<T>
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

      CpmlElectricParam<T> *ptr;
      ptr = static_cast<CpmlElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double by = ptr->b1;
      double bz = ptr->b2;
      double cy = ptr->c1;
      double cz = ptr->c2;
      double kappay = ptr->kappa1;
      double kappaz = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = by * psi1 + cy * (hz(i+1,j+1,k) - hz(i+1,j,k)) / dy;
      psi2 = bz * psi2 + cz * (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;
      
      ex(i,j,k) += dt / eps_inf * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy / kappay -
			       (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz / kappaz +
			       psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  };

  template <typename T> class CpmlEy: public CpmlElectric<T>
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

      CpmlElectricParam<T> *ptr;
      ptr = static_cast<CpmlElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double bz = ptr->b1;
      double bx = ptr->b2;
      double cz = ptr->c1;
      double cx = ptr->c2;
      double kappaz = ptr->kappa1;
      double kappax = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = bz * psi1 + cz * (hx(i,j+1,k+1) - hx(i,j+1,k)) / dz;
      psi2 = bx * psi2 + cx * (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx;
      
      ey(i,j,k) += dt / eps_inf * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz / kappaz -
			       (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx / kappax +
			       psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  };

  template <typename T> class CpmlEz: public CpmlElectric<T>
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

      CpmlElectricParam<T> *ptr;
      ptr = static_cast<CpmlElectricParam<T> *>(parameter);
      double eps_inf = ptr->eps_inf;
      double bx = ptr->b1;
      double by = ptr->b2;
      double cx = ptr->c1;
      double cy = ptr->c2;
      double kappax = ptr->kappa1;
      double kappay = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = bx * psi1 + cx * (hy(i+1,j,k+1) - hy(i,j,k+1)) / dx;
      psi2 = by * psi2 + cy * (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy;
      
      ez(i,j,k) += dt / eps_inf * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx / kappax -
			       (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy / kappay +
			       psi1 - psi2);
    }

  protected:
    using CpmlElectric<T>::param;
  };

  template <typename T> class CpmlMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~CpmlMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<CpmlMagneticParam<T> *>(iter->second);
      }
      param.clear();
    }
    
    PwMaterial<T> *
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::iterator iter = param.find(index);
      if (iter != param.end()) {
      	std::cerr << "Overwriting the existing index." << std::endl;
      	delete static_cast<CpmlMagneticParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const CpmlMagneticParam<T> * const CpmlMagneticParameter_ptr
	= static_cast<const CpmlMagneticParam<T> * const>(parameter);
      CpmlMagneticParam<T> * param_ptr;
      param_ptr = new CpmlMagneticParam<T>();
      param_ptr->mu_inf = CpmlMagneticParameter_ptr->mu_inf;
      param_ptr->b1 = CpmlMagneticParameter_ptr->b1;
      param_ptr->b2 = CpmlMagneticParameter_ptr->b2;
      param_ptr->c1 = CpmlMagneticParameter_ptr->c1;
      param_ptr->c2 = CpmlMagneticParameter_ptr->c2;
      param_ptr->kappa1 = CpmlMagneticParameter_ptr->kappa1;
      param_ptr->kappa2 = CpmlMagneticParameter_ptr->kappa2;
      param_ptr->psi1 = static_cast<T>(0);
      param_ptr->psi2 = static_cast<T>(0);

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class CpmlHx: public CpmlMagnetic<T>
  {
  public:
    void 
    update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   double dy, double dz, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];

      CpmlMagneticParam<T> *ptr;
      ptr = static_cast<CpmlMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double by = ptr->b1;
      double bz = ptr->b2;
      double cy = ptr->c1;
      double cz = ptr->c2;
      double kappay = ptr->kappa1;
      double kappaz = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = by * psi1 + cy * (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy;
      psi2 = bz * psi2 + cz * (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz;

      hx(i,j,k) -= dt / mu_inf * ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy / kappay -
			      (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz / kappaz +
			      psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::param;
  };

  template <typename T> class CpmlHy: public CpmlMagnetic<T>
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

      CpmlMagneticParam<T> *ptr;
      ptr = static_cast<CpmlMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double bz = ptr->b1;
      double bx = ptr->b2;
      double cz = ptr->c1;
      double cx = ptr->c2;
      double kappaz = ptr->kappa1;
      double kappax = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = bz * psi1 + cz * (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz;
      psi2 = bx * psi2 + cx * (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx;

      hy(i,j,k) -= dt / mu_inf * ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz / kappaz -
			      (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx / kappax +
			      psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::param;
  };

  template <typename T> class CpmlHz: public CpmlMagnetic<T>
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

      CpmlMagneticParam<T> *ptr;
      ptr = static_cast<CpmlMagneticParam<T> *>(parameter);
      double mu_inf = ptr->mu_inf;
      double bx = ptr->b1;
      double by = ptr->b2;
      double cx = ptr->c1;
      double cy = ptr->c2;
      double kappax = ptr->kappa1;
      double kappay = ptr->kappa2;
      T& psi1 = ptr->psi1;
      T& psi2 = ptr->psi2;

      psi1 = bx * psi1 + cx * (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx;
      psi2 = by * psi2 + cy * (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy;
      
      hz(i,j,k) -= dt / mu_inf * ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx / kappax -
				  (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy / kappay +
				  psi1 - psi2);
    }

  protected:
    using CpmlMagnetic<T>::param;
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif // PW_CPML_HH_
