#ifndef PW_DUMMY_HH_
#define PW_DUMMY_HH_

#include <iostream>
#include "pw_material.hh"

namespace gmes
{
  template <typename T> struct DummyElectricParam: public ElectricParam<T>
  {
  };

  template <typename T> struct DummyMagneticParam: public MagneticParam<T>
  {
  };
  
  template <typename T> class DummyElectric: public MaterialElectric<T>
  {
  public:
    ~DummyElectric()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DummyElectricParam<T> *>(iter->second);
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
	delete static_cast<DummyElectricParam<T> *>(iter->second);
	param.erase(iter);
      }

      DummyElectricParam<T> *param_ptr;
      param_ptr = new DummyElectricParam<T>();
      param_ptr->eps_inf = static_cast<const DummyElectricParam<T> *>(parameter)->eps_inf;
	  
      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

    void 
    update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
    }

  protected:
    using MaterialElectric<T>::param;
  };
  
  template <typename T> class DummyEx: public DummyElectric<T>
  {
  };

  template <typename T> class DummyEy: public DummyElectric<T>
  {
  };

  template <typename T> class DummyEz: public DummyElectric<T>
  {
  };

  template <typename T> class DummyMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~DummyMagnetic()
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	delete static_cast<DummyMagneticParam<T> *>(iter->second);
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
	delete static_cast<DummyMagneticParam<T> *>(iter->second);
	param.erase(iter);
      }
      
      DummyMagneticParam<T> *param_ptr;
      param_ptr = new DummyMagneticParam<T>();
      param_ptr->mu_inf = static_cast<const DummyMagneticParam<T> *>(parameter)->mu_inf;
      
      param.insert(std::make_pair(index, param_ptr));

      return this;
    }

    void 
    update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n,
	   const int* const idx, int idx_size, 
	   PwMaterialParam * const parameter)
    {
    }

  protected:
    using MaterialMagnetic<T>::param;
  };

  template <typename T> class DummyHx: public DummyMagnetic<T>
  {
  };

  template <typename T> class DummyHy: public DummyMagnetic<T>
  {
  };

  template <typename T> class DummyHz: public DummyMagnetic<T>
  {
  };
}

#endif /*PW_DUMMY_HH_*/
