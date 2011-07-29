#ifndef PW_CONST_HH_
#define PW_CONST_HH_

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
      for(MapType::const_iterator iter = param.begin(); iter != param.end(); iter++) {
	delete[] iter->first;
	delete static_cast<ConstElectricParam *>(iter->second);
	}
      param.clear();
    }
    
    void 
    attach(const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      int *idx_ptr = new int[3];
      std::copy(idx, idx + idx_size, idx_ptr);

      ConstElectricParam *param_ptr;
      param_ptr = new ConstElectricParam();
      param_ptr->eps = static_cast<ConstElectricParam *>(parameter)->eps;
      param_ptr->value = static_cast<ConstMagneticParam *>(parameter)->value;
      param.insert(std::make_pair(idx_ptr, param_ptr));
    }
    
    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double value = static_cast<ConstElectricParam *>(parameter)->value;

      inplace_field(i,j,k) = value;
    }

  protected:
    using PwMaterial<T>::param;
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
      for(MapType::const_iterator iter = param.begin(); iter != param.end(); iter++) {
	delete[] iter->first;
	delete static_cast<ConstMagneticParam *>(iter->second);
	}
      param.clear();
    }
    
    void 
    attach(const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      int *idx_ptr = new int[3];
      std::copy(idx, idx + idx_size, idx_ptr);

      ConstMagneticParam *param_ptr;
      param_ptr = new ConstMagneticParam();
      param_ptr->mu = static_cast<ConstMagneticParam *>(parameter)->mu;
      param_ptr->value = static_cast<ConstMagneticParam *>(parameter)->value;

      param.insert(std::make_pair(idx_ptr, param_ptr));
    }

    void 
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int idx[3], int idx_size, 
	   const PwMaterialParam * const parameter)
    {
      int i = idx[0], j = idx[1], k = idx[2];
      double value = static_cast<ConstElectricParam *>(parameter)->value;

      inplace_field(i,j,k) = value;
    }

  protected:
    using PwMaterial<T>::param;
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
