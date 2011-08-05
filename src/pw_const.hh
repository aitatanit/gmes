#ifndef PW_CONST_HH_
#define PW_CONST_HH_

#include <iostream>
#include <utility>
#include "pw_material.hh"

#define inplace_field(i,j,k) inplace_field[((i)*inplace_dim2+(j))*inplace_dim3+(k)]

namespace gmes
{
  struct ConstElectricParam: public ElectricParam 
  {
    double value;
  };
    
  struct ConstMagneticParam: public MagneticParam 
  {
    double value;
  };

  template <typename T> class ConstElectric: public MaterialElectric<T>
  {
  public:
    ~ConstElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<ConstElectricParam *>(iter->second);
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
	delete static_cast<ConstElectricParam *>(iter->second);
	param.erase(iter);
      }

      const ConstElectricParam* ConstElectricParameter_ptr
	= static_cast<const ConstElectricParam *>(parameter);      
      ConstElectricParam *param_ptr = new ConstElectricParam();
      param_ptr->eps = ConstElectricParameter_ptr->eps;
      param_ptr->value = ConstElectricParameter_ptr->value;

      param.insert(std::make_pair(index, param_ptr));
    }
    
    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double value = static_cast<const ConstElectricParam *>(parameter)->value;

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
	delete static_cast<ConstMagneticParam *>(iter->second);
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
	delete static_cast<ConstMagneticParam *>(iter->second);
	param.erase(iter);
      }

      const ConstMagneticParam* ConstMagneticParameter_ptr
	= static_cast<const ConstMagneticParam *>(parameter);
      ConstMagneticParam *param_ptr = new ConstMagneticParam();
      param_ptr->mu = ConstMagneticParameter_ptr->mu;
      param_ptr->value = ConstMagneticParameter_ptr->value;

      param.insert(std::make_pair(index, param_ptr));
    }

    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double value = static_cast<const ConstElectricParam *>(parameter)->value;

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
