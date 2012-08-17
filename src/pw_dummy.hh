#ifndef PW_DUMMY_HH_
#define PW_DUMMY_HH_

#include <iostream>
#include "pw_material.hh"

namespace gmes
{
  template <typename T> struct DummyElectricParam: public ElectricParam<T>
  {
  }; // template DummyElectricParam

  template <typename T> struct DummyMagneticParam: public MagneticParam<T>
  {
  }; // template DummyMagneticParam
  
  template <typename T> class DummyElectric: public MaterialElectric<T>
  {
  public:
    ~DummyElectric()
    {
      for (auto v: param) {
	delete static_cast<DummyElectricParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T> *
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto iter = param.find(index);
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
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      // Dummy does nothing.
    }

  protected:
    using MaterialElectric<T>::param;
  }; // template DummyElectric
  
  template <typename T> class DummyEx: public DummyElectric<T>
  {
  }; // template DummyEx

  template <typename T> class DummyEy: public DummyElectric<T>
  {
  }; // template DummyEy

  template <typename T> class DummyEz: public DummyElectric<T>
  {
  }; // template DummyEz

  template <typename T> class DummyMagnetic: public MaterialMagnetic<T>
  {
  public:
    ~DummyMagnetic()
    {
      for (auto v: param) {
	delete static_cast<DummyMagneticParam<T> *>(v.second);
      }
      param.clear();
    }
    
    PwMaterial<T> *
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam * const parameter)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto iter = param.find(index);
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
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      // Dummy does nothing.
    }

  protected:
    using MaterialMagnetic<T>::param;
  }; // template DummyMagnetic

  template <typename T> class DummyHx: public DummyMagnetic<T>
  {
  }; // template DummyHx

  template <typename T> class DummyHy: public DummyMagnetic<T>
  {
  }; // template DummyHy

  template <typename T> class DummyHz: public DummyMagnetic<T>
  {
  }; // template DummyHz
} // namespace gmes

#endif // PW_DUMMY_HH_
