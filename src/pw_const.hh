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
  };
    
  template <typename T> struct ConstMagneticParam: public MagneticParam<T>
  {
    T value;
  };

  template <typename T> class ConstElectric: public MaterialElectric<T>
  {
  public:
    ~ConstElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<ConstElectricParam<T> *>(iter->second);
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
      	delete static_cast<ConstElectricParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const ConstElectricParam<T> *ConstElectricParameter_ptr
	= static_cast<const ConstElectricParam<T> *>(parameter);      
      ConstElectricParam<T> *param_ptr = new ConstElectricParam<T>();

      param_ptr->eps_inf = ConstElectricParameter_ptr->eps_inf;
      param_ptr->value = ConstElectricParameter_ptr->value;

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }
    
    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const int* const idx, int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      const T& value 
	= static_cast<const ConstElectricParam<T> *>(parameter)->value;

      inplace_field(i,j,k) = value;
    }

  protected:
    using MaterialElectric<T>::param;
  };

  template <typename T> class ConstEx: public ConstElectric<T>
  {
  };

  template <typename T> class ConstEy: public ConstElectric<T>
  {
  };

  template <typename T> class ConstEz: public ConstElectric<T>
  {
  };

  template <typename T> class ConstMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~ConstMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<ConstMagneticParam<T> *>(iter->second);
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
      	delete static_cast<ConstMagneticParam<T> *>(iter->second);
      	param.erase(iter);
      }

      const ConstMagneticParam<T> *ConstMagneticParameter_ptr
	= static_cast<const ConstMagneticParam<T> *>(parameter);
      ConstMagneticParam<T> *param_ptr = new ConstMagneticParam<T>();

      param_ptr->mu_inf = ConstMagneticParameter_ptr->mu_inf;
      param_ptr->value = ConstMagneticParameter_ptr->value;

      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int* const idx, int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      const T& value 
	= static_cast<const ConstElectricParam<T> *>(parameter)->value;

      inplace_field(i,j,k) = value;
    }

  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class ConstHx: public ConstMagnetic<T>
  {
  };

  template <typename T> class ConstHy: public ConstMagnetic<T>
  {
  };

  template <typename T> class ConstHz: public ConstMagnetic<T>
  {
  };
}

#undef inplace_field

#endif /*PW_CONST_HH_*/
