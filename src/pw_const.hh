#ifndef PW_CONST_HH_
#define PW_CONST_HH_

#include <iostream>
#include <utility>
#include "pw_material.hh"

#define inplace_field(i,j,k) inplace_field[inplace_dim1==1&&inplace_dim2==1&&inplace_dim3==1?0:((i)*inplace_dim2+(j))*inplace_dim3+(k)]

namespace gmes
{
  template <typename T> struct ConstElectricParam: public ElectricParam<T>
  {
    T value;
  }; // template ConstElectricParam
    
  template <typename T> struct ConstMagneticParam: public MagneticParam<T>
  {
    T value;
  }; // template ConstMagneticParam

  template <typename T> class ConstElectric: public MaterialElectric<T>
  {
  public:
    ~ConstElectric()
    {
      for (auto v: param) {
	delete static_cast<ConstElectricParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto const_param_ptr = static_cast<const ConstElectricParam<T>*>(pm_param_ptr);
      auto new_param_ptr = new ConstElectricParam<T>();

      new_param_ptr->eps_inf = const_param_ptr->eps_inf;
      new_param_ptr->value = const_param_ptr->value;

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }
    
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
	const auto const_param_ptr = static_cast<ConstElectricParam<T>*>(v.second);
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, *const_param_ptr);
      }
    }

  private:
    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const Index3& idx, 
	   const ConstElectricParam<T>& const_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      inplace_field(i,j,k) = const_param.value;
    }

  protected:
    using MaterialElectric<T>::param;
  }; // template ConstElectric

  template <typename T> class ConstEx: public ConstElectric<T>
  {
  }; // template ConstEx

  template <typename T> class ConstEy: public ConstElectric<T>
  {
  }; // template ConstEy

  template <typename T> class ConstEz: public ConstElectric<T>
  {
  }; // template ConstEz

  template <typename T> class ConstMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~ConstMagnetic()
    {
      for (auto v: param) {
	delete static_cast<ConstMagneticParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T>*
    attach(const int* const idx, int idx_size, 
	   const PwMaterialParam* const pm_param_ptr)
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto const_param_ptr = static_cast<const ConstMagneticParam<T>*>(pm_param_ptr);
      auto new_param_ptr = new ConstMagneticParam<T>();

      new_param_ptr->mu_inf = const_param_ptr->mu_inf;
      new_param_ptr->value = const_param_ptr->value;

      param.insert(std::make_pair(index, new_param_ptr));

      return this;
    }

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
	const auto const_param_ptr = static_cast<ConstMagneticParam<T>*>(v.second);
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, 
    	       v.first, *const_param_ptr);
      }
    }

  private:
    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const Index3& idx, 
	   const ConstMagneticParam<T>& const_param) const
    {
      const int i = idx[0], j = idx[1], k = idx[2];
      inplace_field(i,j,k) = const_param.value;
    }

  protected:
    using MaterialMagnetic<T>::param;
  }; // template ConstMagnetic

  template <typename T> class ConstHx: public ConstMagnetic<T>
  {
  }; // template ConstHx

  template <typename T> class ConstHy: public ConstMagnetic<T>
  {
  }; // template ConstHy

  template <typename T> class ConstHz: public ConstMagnetic<T>
  {
  }; // template ConstHz
} // namespace gmes

#undef inplace_field

#endif // PW_CONST_HH_
