#ifndef PW_CONST_HH_
#define PW_CONST_HH_

#include <iostream>
#include <utility>
#include "pw_material.hh"

#define inplace_field(i,j,k) inplace_field[inplace_dim1==1&&inplace_dim2==1&&inplace_dim3==1?0:((i)*inplace_dim2+(j))*inplace_dim3+(k)]

namespace gmes
{
  template <typename T> 
  struct ConstElectricParam: public ElectricParam<T>
  {
    T value;
  }; // template ConstElectricParam
    
  template <typename T> 
  struct ConstMagneticParam: public MagneticParam<T>
  {
    T value;
  }; // template ConstMagneticParam

  template <typename T> 
  class ConstElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return ConstElectric<T>::tag;
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

      const auto& const_param = *static_cast<const ConstElectricParam<T>*>(pm_param_ptr);
      
      idx_list.push_back(index);
      param_list.push_back(const_param);

      return this;
    }
    
    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto const_ptr = static_cast<const ConstElectric<T>*>(pm_ptr);
      std::copy(const_ptr->idx_list.begin(), const_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(const_ptr->param_list.begin(), const_ptr->param_list.end(), std::back_inserter(param_list));
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
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, *idx, *param);
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
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<ConstElectricParam<T> > param_list;

  private:
    static const std::string tag; // "ConstElectric"
  }; // template ConstElectric

  template <typename T>
  const std::string ConstElectric<T>::tag = "ConstElectric";

  template <typename T> 
  class ConstEx: public ConstElectric<T>
  {
  }; // template ConstEx

  template <typename T> 
  class ConstEy: public ConstElectric<T>
  {
  }; // template ConstEy

  template <typename T> 
  class ConstEz: public ConstElectric<T>
  {
  }; // template ConstEz

  template <typename T> 
  class ConstMagnetic: public MaterialMagnetic<T>
  {
  public:
    const std::string& 
    name() const
    {
      return ConstMagnetic<T>::tag;
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

      const auto& const_param = *static_cast<const ConstMagneticParam<T>*>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(const_param);

      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto const_ptr = static_cast<const ConstMagnetic<T>*>(pm_ptr);
      std::copy(const_ptr->idx_list.begin(), const_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(const_ptr->param_list.begin(), const_ptr->param_list.end(), std::back_inserter(param_list));
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
      for (auto idx = idx_list.begin(), param = param_list.begin();
	   idx != idx_list.end(); ++idx, ++param) {
    	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
    	       in_field1, in1_dim1, in1_dim2, in1_dim3,
    	       in_field2, in2_dim1, in2_dim2, in2_dim3,
    	       d1, d2, dt, n, *idx, *param);
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
    using MaterialMagnetic<T>::position;
    using MaterialMagnetic<T>::idx_list;
    std::vector<ConstMagneticParam<T> > param_list;

  private:
    static const std::string tag; // "ConstMagnetic"
  }; // template ConstMagnetic

  template <typename T>
  const std::string ConstMagnetic<T>::tag = "ConstMagnetic";

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
