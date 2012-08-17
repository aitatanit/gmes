/* The implementation of Drude-critical points model based on the following
 * articles.
 *
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "An analytic model for the 
 * optical properties of gold," The Journal of Chemical Physics, vol. 125,
 * no. 16, pp. 164705-3, Oct. 2006.
 *
 * P. G. Etchegoin, E. C. Le Ru, and M. Meyer, "Erratum: An analytic model
 * for the optical properties of gold" [J. Chem. Phys. 125, 164705 (2006)]," 
 * The Journal of Chemical Physics, vol. 127, no. 18, pp. 189901-1, Nov. 2007.
 *
 * A. Taflove and S. C. Hagness, Computational Electrodynamics: The Finite-
 * Difference Time-Domain Method, 3rd ed. Artech House Publishers, 2005.
 */

#ifndef PW_DCP_HH_
#define PW_DCP_HH_

#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include "pw_dielectric.hh"

#define ex(i,j,k) ex[ex_y_size==1?0:((i)*ex_y_size+(j))*ex_z_size+(k)]
#define ey(i,j,k) ey[ey_z_size==1?0:((i)*ey_y_size+(j))*ey_z_size+(k)]
#define ez(i,j,k) ez[ez_x_size==1?0:((i)*ez_y_size+(j))*ez_z_size+(k)]
#define hx(i,j,k) hx[hx_y_size==1?0:((i)*hx_y_size+(j))*hx_z_size+(k)]
#define hy(i,j,k) hy[hy_z_size==1?0:((i)*hy_y_size+(j))*hy_z_size+(k)]
#define hz(i,j,k) hz[hz_x_size==1?0:((i)*hz_y_size+(j))*hz_z_size+(k)]

// The classes should be rewritten using template specialization
// to increase the calculation speed.
namespace gmes
{
  /**************************************************/
  /* Auxiliary Differential Equation Implementation */
  /**************************************************/

  template <typename T> struct DcpAdeElectricParam: ElectricParam<T>
  {
    // parameters for the ADE of the Drude model
    std::vector<std::array<double, 3> > a;
    // parameters for the ADE of critical points model
    std::vector<std::array<double, 4> > b;
    // parameters for the electric field update equations
    std::array<double, 4> c;
    T e_old;
    std::vector<T> q_old, q_now, p_old, p_now;
  };
  
  template <typename T> struct DcpAdeMagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DcpAdeElectric: public MaterialElectric<T>
  {
  public:
    ~DcpAdeElectric()
    {
      for (auto v: param) {
	delete static_cast<DcpAdeElectricParam<T> *>(v.second);
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
      	delete static_cast<DcpAdeElectricParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const DcpAdeElectricParam<T> * const DcpAdeElectricParameter_ptr
	= static_cast<const DcpAdeElectricParam<T> * const>(parameter);
      DcpAdeElectricParam<T> *param_ptr;
      param_ptr = new DcpAdeElectricParam<T>();

      param_ptr->eps_inf = DcpAdeElectricParameter_ptr->eps_inf;
      std::copy(DcpAdeElectricParameter_ptr->a.begin(),
		DcpAdeElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(DcpAdeElectricParameter_ptr->b.begin(),
		DcpAdeElectricParameter_ptr->b.end(),
		std::back_inserter(param_ptr->b));
      std::copy(DcpAdeElectricParameter_ptr->c.begin(),
		DcpAdeElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->e_old = static_cast<T>(0);
      param_ptr->q_old.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->q_now.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->p_old.resize(param_ptr->b.size(), static_cast<T>(0));
      param_ptr->p_now.resize(param_ptr->b.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }
    
    T 
    dps_sum(const T& init, const DcpAdeElectricParam<T> * const ptr) const
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      const std::vector<T>& q_old = ptr->q_old;
      const std::vector<T>& q_now = ptr->q_now;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	sum += a[i][0] * q_old[i] + (a[i][1] - 1) * q_now[i];
      }

      return sum;
    }
      
    T 
    cps_sum(const T& init, const DcpAdeElectricParam<T> * const ptr) const
    {
      const std::vector<std::array<double, 4> >& b = ptr->b;
      const std::vector<T>& p_old = ptr->p_old;
      const std::vector<T>& p_now = ptr->p_now;

      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	sum += b[i][0] * p_old[i] + (b[i][1] - 1) * p_now[i];
      }
      
      return sum;
    }

    void 
    update_q(const T& e_now, DcpAdeElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      std::vector<T>& q_old = ptr->q_old;
      std::vector<T>& q_now = ptr->q_now;

      std::vector<T> q_new(a.size());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * e_now;
      }
      
      std::copy(q_now.begin(), q_now.end(), q_old.begin());
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
    }
    
    void 
    update_p(const T& e_old, const T& e_now, const T& e_new,
	     DcpAdeElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 4> >& b = ptr->b;
      std::vector<T>& p_old = ptr->p_old;
      std::vector<T>& p_now = ptr->p_now;
    
      std::vector<T> p_new(b.size());
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] 
	  * (e_old - e_new) + b[i][3] * e_now;
      }
      std::copy(p_now.begin(), p_now.end(), p_old.begin());
      std::copy(p_new.begin(), p_new.end(), p_now.begin());
    }
  
  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class DcpAdeEx: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;
 
      const T& e_now = ex(i,j,k);
      const T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			      (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) 
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ex(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeEy: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;
      
      const T& e_now = ey(i,j,k);
      T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;

      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ey(i,j,k) = e_new;
    }
    
  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeEz: public DcpAdeElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAdeElectricParam<T> *ptr;
      ptr = static_cast<DcpAdeElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;

      const T& e_now = ez(i,j,k);
      T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_now, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ez(i,j,k) = e_new;
    }

  protected:
    using DcpAdeElectric<T>::param;
    using DcpAdeElectric<T>::dps_sum;
    using DcpAdeElectric<T>::cps_sum;
    using DcpAdeElectric<T>::update_q;
    using DcpAdeElectric<T>::update_p;
  };

  template <typename T> class DcpAdeHx: public DielectricHx<T>
  {
  };

  template <typename T> class DcpAdeHy: public DielectricHy<T>
  {
  };

  template <typename T> class DcpAdeHz: public DielectricHz<T>
  {
  };

  /***********************************************************/
  /* Improved Auxiliary Differential Equation Implementation */
  /***********************************************************/

  typedef std::vector<std::array<double, 3> > Ade2CoeffA;
  typedef std::vector<std::array<double, 5> > Ade2CoeffB;
  typedef std::array<double, 4> Ade2CoeffC;
    
  template <typename T> struct DcpAde2ElectricParam: ElectricParam<T>
  {
    // parameters for the ADE of the Drude terms
    Ade2CoeffA a;
    // parameters for the ADE of critical point terms
    Ade2CoeffB b;
    // parameters for the electric field update equations
    Ade2CoeffC c;
    T e_old;
    std::vector<T> q_old, q_now, p_old, p_now;
  };
  
  template <typename T> struct DcpAde2MagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DcpAde2Electric: public MaterialElectric<T>
  {
  public:
    ~DcpAde2Electric()
    {
      for (auto v: param) {
	delete static_cast<DcpAde2ElectricParam<T> *>(v.second);
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
      	delete static_cast<DcpAde2ElectricParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const DcpAde2ElectricParam<T> * const DcpAde2ElectricParameter_ptr
	= static_cast<const DcpAde2ElectricParam<T> * const>(parameter);
      DcpAde2ElectricParam<T> *param_ptr;
      param_ptr = new DcpAde2ElectricParam<T>();

      param_ptr->eps_inf = DcpAde2ElectricParameter_ptr->eps_inf;
      std::copy(DcpAde2ElectricParameter_ptr->a.begin(),
		DcpAde2ElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(DcpAde2ElectricParameter_ptr->b.begin(),
		DcpAde2ElectricParameter_ptr->b.end(),
		std::back_inserter(param_ptr->b));
      std::copy(DcpAde2ElectricParameter_ptr->c.begin(),
		DcpAde2ElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->e_old = static_cast<T>(0);
      param_ptr->q_old.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->q_now.resize(param_ptr->a.size(), static_cast<T>(0));
      param_ptr->p_old.resize(param_ptr->b.size(), static_cast<T>(0));
      param_ptr->p_now.resize(param_ptr->b.size(), static_cast<T>(0));

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }
    
    T 
    dps_sum(const T& init, const DcpAde2ElectricParam<T> * const ptr) const
    {
      const Ade2CoeffA& a = ptr->a;
      const std::vector<T>& q_old = ptr->q_old;
      const std::vector<T>& q_now = ptr->q_now;
      
      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	sum += (1 - a[i][1]) * q_now[i] - a[i][0] * q_old[i];
      }

      return sum;
    }
    
    T 
    cps_sum(const T& init, const DcpAde2ElectricParam<T> * const ptr) const
    {
      const Ade2CoeffB& b = ptr->b;
      const std::vector<T>& p_old = ptr->p_old;
      const std::vector<T>& p_now = ptr->p_now;

      T sum(init);
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	sum += (1 - b[i][1]) * p_now[i] - b[i][0] * p_old[i];
      }
      
      return sum;
    }

    void 
    update_q(const T& e_old, const T& e_now, const T& e_new,
	     DcpAde2ElectricParam<T> * const ptr)
    {
      const Ade2CoeffA& a = ptr->a;
      std::vector<T>& q_old = ptr->q_old;
      std::vector<T>& q_now = ptr->q_now;

      std::vector<T> q_new(a.size());
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	q_new[i] = a[i][0] * q_old[i] + a[i][1] * q_now[i] + a[i][2] * (e_old + 2.0 * e_now + e_new);
      }
      
      std::copy(q_now.begin(), q_now.end(), q_old.begin());
      std::copy(q_new.begin(), q_new.end(), q_now.begin());
    }
    
    void 
    update_p(const T& e_old, const T& e_now, const T& e_new,
	     DcpAde2ElectricParam<T> * const ptr)
    {
      const Ade2CoeffB& b = ptr->b;
      std::vector<T>& p_old = ptr->p_old;
      std::vector<T>& p_now = ptr->p_now;
    
      std::vector<T> p_new(b.size());
      for (typename std::vector<T>::size_type i = 0; i < b.size(); ++i) {
	p_new[i] = b[i][0] * p_old[i] + b[i][1] * p_now[i] + b[i][2] 
	  * e_old + b[i][3] * e_now + b[i][4] * e_new;
      }
      std::copy(p_now.begin(), p_now.end(), p_old.begin());
      std::copy(p_new.begin(), p_new.end(), p_now.begin());
    }
  
  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class DcpAde2Ex: public DcpAde2Electric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAde2ElectricParam<T> *ptr;
      ptr = static_cast<DcpAde2ElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;

      const T& e_now = ex(i,j,k);
      const T e_new = c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
			      (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) 
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_old, e_now, e_new, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ex(i,j,k) = e_new;
    }
    
  protected:
    using DcpAde2Electric<T>::param;
    using DcpAde2Electric<T>::dps_sum;
    using DcpAde2Electric<T>::cps_sum;
    using DcpAde2Electric<T>::update_q;
    using DcpAde2Electric<T>::update_p;
  };

  template <typename T> class DcpAde2Ey: public DcpAde2Electric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAde2ElectricParam<T> *ptr;
      ptr = static_cast<DcpAde2ElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;
      
      const T& e_now = ey(i,j,k);
      T e_new = c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
			(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;

      update_q(e_old, e_now, e_new, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ey(i,j,k) = e_new;
    }
    
  protected:
    using DcpAde2Electric<T>::param;
    using DcpAde2Electric<T>::dps_sum;
    using DcpAde2Electric<T>::cps_sum;
    using DcpAde2Electric<T>::update_q;
    using DcpAde2Electric<T>::update_p;
  };

  template <typename T> class DcpAde2Ez: public DcpAde2Electric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpAde2ElectricParam<T> *ptr;
      ptr = static_cast<DcpAde2ElectricParam<T> *>(parameter);
      const std::array<double, 4>& c = ptr->c;
      T& e_old = ptr->e_old;

      const T& e_now = ez(i,j,k);
      T e_new = c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
			(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy)
	+ c[1] * (dps_sum(static_cast<T>(0), ptr) + 
		  cps_sum(static_cast<T>(0), ptr))
	+ c[2] * e_old + c[3] * e_now;
      
      update_q(e_old, e_now, e_new, ptr);
      update_p(e_old, e_now, e_new, ptr);
      
      e_old = e_now;
      ez(i,j,k) = e_new;
    }

  protected:
    using DcpAde2Electric<T>::param;
    using DcpAde2Electric<T>::dps_sum;
    using DcpAde2Electric<T>::cps_sum;
    using DcpAde2Electric<T>::update_q;
    using DcpAde2Electric<T>::update_p;
  };

  template <typename T> class DcpAde2Hx: public DielectricHx<T>
  {
  };

  template <typename T> class DcpAde2Hy: public DielectricHy<T>
  {
  };

  template <typename T> class DcpAde2Hz: public DielectricHz<T>
  {
  };

  /*********************************************************/
  /* Piecewise-Linear Recursive Convolution Implementation */
  /*********************************************************/

  template <typename T> struct DcpPlrcElectricParam: ElectricParam<T>
  {
    std::vector<std::array<double, 3> > a;
    std::vector<std::array<std::complex<double>, 3> > b;
    std::array<double, 3> c;
    // *_re and *_im are for the real and imaginary part of the e-field, 
    // respectively.
    std::vector<double> psi_dp_re, psi_dp_im;
    std::vector<std::complex<double> > psi_cp_re, psi_cp_im;
  };

  template <typename T> struct DcpPlrcMagneticParam: MagneticParam<T>
  {
  };

  template <typename T> class DcpPlrcElectric: public MaterialElectric<T>
  {
  public:
    ~DcpPlrcElectric()
    {
      for (auto v: param) {
	delete static_cast<DcpPlrcElectricParam<T> *>(v.second);
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
      	delete static_cast<DcpPlrcElectricParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const DcpPlrcElectricParam<T> * const DcpPlrcElectricParameter_ptr
	= static_cast<const DcpPlrcElectricParam<T> * const>(parameter);
      DcpPlrcElectricParam<T> *param_ptr;
      param_ptr = new DcpPlrcElectricParam<T>();

      param_ptr->eps_inf = DcpPlrcElectricParameter_ptr->eps_inf;
      std::copy(DcpPlrcElectricParameter_ptr->a.begin(),
		DcpPlrcElectricParameter_ptr->a.end(),
		std::back_inserter(param_ptr->a));
      std::copy(DcpPlrcElectricParameter_ptr->b.begin(),
		DcpPlrcElectricParameter_ptr->b.end(),
		std::back_inserter(param_ptr->b));
      std::copy(DcpPlrcElectricParameter_ptr->c.begin(),
		DcpPlrcElectricParameter_ptr->c.end(),
		param_ptr->c.begin());
      param_ptr->psi_dp_re.resize(param_ptr->a.size(), 0);
      param_ptr->psi_dp_im.resize(param_ptr->a.size(), 0);
      param_ptr->psi_cp_re.resize(param_ptr->b.size(), std::complex<double>(0));
      param_ptr->psi_cp_im.resize(param_ptr->b.size(), std::complex<double>(0));
      
      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

    void 
    update_psi_dp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T> * const ptr)
    {
      const std::vector<std::array<double, 3> >& a = ptr->a;
      std::vector<double>& psi_dp_re = ptr->psi_dp_re;
      std::vector<double>& psi_dp_im = ptr->psi_dp_im;
      
      for (typename std::vector<T>::size_type i = 0; i < a.size(); ++i) {
	psi_dp_re[i] = a[i][0] * e_new.real() + a[i][1] * e_now.real() 
	  + a[i][2] * psi_dp_re[i];
	psi_dp_im[i] = a[i][0] * e_new.imag() + a[i][1] * e_now.imag() 
	  + a[i][2] * psi_dp_im[i];
      }
    }

    void 
    update_psi_cp(const std::complex<double>& e_now, 
		  const std::complex<double>& e_new,
		  DcpPlrcElectricParam<T> * const ptr)
    {
      const std::vector<std::array<std::complex<double>, 3> >& b = ptr->b;
      std::vector<std::complex<double> >& psi_cp_re = ptr->psi_cp_re;
      std::vector<std::complex<double> >& psi_cp_im = ptr->psi_cp_im;
      
      for (typename std::vector<std::complex<double> >::size_type i = 0; 
	   i < b.size(); ++i) {
	psi_cp_re[i] = b[i][0] * e_new.real() + b[i][1] * e_now.real()
	  + b[i][2] * psi_cp_re[i];
	psi_cp_im[i] = b[i][0] * e_new.imag() + b[i][1] * e_now.imag()
	  + b[i][2] * psi_cp_im[i];
      }
    }

    std::complex<double> 
    psi_total(const DcpPlrcElectricParam<T> * const ptr) const
    {
      const std::vector<double>& psi_dp_re = ptr->psi_dp_re;
      const std::vector<double>& psi_dp_im = ptr->psi_dp_im;
      const std::vector<std::complex<double> >& psi_cp_re = ptr->psi_cp_re;
      const std::vector<std::complex<double> >& psi_cp_im = ptr->psi_cp_im;

      double psi_re = 
	std::accumulate(psi_dp_re.begin(), psi_dp_re.end(), 0.0) + 
	std::accumulate(psi_cp_re.begin(), psi_cp_re.end(), std::complex<double>(0)).real();
      
      double psi_im = 
	std::accumulate(psi_dp_im.begin(), psi_dp_im.end(), 0.0) +
	std::accumulate(psi_cp_im.begin(), psi_cp_im.end(), std::complex<double>(0)).real();
      
      return std::complex<double>(psi_re, psi_im);
    }
    
  protected:
    using MaterialElectric<T>::param;
  };
  
  template <typename S, typename T>
  static inline T& 
  assign(const std::complex<S>& in, T& out)
  {
    return out = static_cast<T>(in.real());
  }

  template <typename S, typename T>
  static inline std::complex<T>& 
  assign(const std::complex<S>& in, std::complex<T>& out)
  {
    return out = in;
  }

  template <typename T> class DcpPlrcEx: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T* const ex, int ex_x_size, int ex_y_size, int ex_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   double dy, double dz, double dt, double n,
	   const Index3& idx, 
	   PwMaterialParam* const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ex(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy - 
		(hy(i+1,j,k+1) - hy(i+1,j,k)) / dz) +
	c[1] * e_now + c[2] * psi_total(ptr);
      
      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);

      assign(e_new, ex(i,j,k));
  }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcEy: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
	   double dz, double dx, double dt, double n,
	   const Index3& idx, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ey(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz - 
		(hz(i+1,j+1,k) - hz(i,j+1,k)) / dx) +
	c[1] * e_now + c[2] * psi_total(ptr);

      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);

      assign(e_new, ey(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcEz: public DcpPlrcElectric<T>
  {
  public:
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for (auto v: param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, v.second);
      }
    }

  private:
    void 
    update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
	   const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
	   const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
	   double dx, double dy, double dt, double n,
	   const Index3& idx, 
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      
      DcpPlrcElectricParam<T> *ptr;
      ptr = static_cast<DcpPlrcElectricParam<T> *>(parameter);
      const std::array<double, 3>& c = ptr->c;

      const std::complex<double> e_now = ez(i,j,k);
      const std::complex<double> e_new = 
	c[0] * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx - 
		(hx(i,j+1,k+1) - hx(i,j,k+1)) / dy) +
	c[1] * e_now + c[2] * psi_total(ptr);

      update_psi_dp(e_now, e_new, ptr);
      update_psi_cp(e_now, e_new, ptr);
      
      assign(e_new, ez(i,j,k));
    }

  protected:
    using DcpPlrcElectric<T>::param;
    using DcpPlrcElectric<T>::update_psi_dp;
    using DcpPlrcElectric<T>::update_psi_cp;
    using DcpPlrcElectric<T>::psi_total;
  };

  template <typename T> class DcpPlrcHx: public DielectricHx<T>
  {
  };

  template <typename T> class DcpPlrcHy: public DielectricHy<T>
  {
  };

  template <typename T> class DcpPlrcHz: public DielectricHz<T>
  {
  };
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DCP_HH_*/
