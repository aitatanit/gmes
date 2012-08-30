#ifndef PW_DUMMY_HH_
#define PW_DUMMY_HH_

#include "pw_material.hh"

namespace gmes
{
  template <typename T> 
  struct DummyElectricParam: public ElectricParam<T>
  {
  }; // template DummyElectricParam

  template <typename T> 
  struct DummyMagneticParam: public MagneticParam<T>
  {
  }; // template DummyMagneticParam
  
  template <typename T> 
  class DummyElectric: public MaterialElectric<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DummyElectric<T>::tag;
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
      
      const auto& dummy_param = *static_cast<const DummyElectricParam<T>*>(pm_param_ptr);

      idx_list.push_back(index);
      param_list.push_back(dummy_param);

      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dummy_ptr = static_cast<const DummyElectric<T>*>(pm_ptr);
      std::copy(dummy_ptr->idx_list.begin(), dummy_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(dummy_ptr->param_list.begin(), dummy_ptr->param_list.end(), std::back_inserter(param_list));
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
    using MaterialElectric<T>::position;
    using MaterialElectric<T>::idx_list;
    std::vector<DummyElectricParam<T> > param_list;

  private:
    static const std::string tag; // "DrudeElectric"
  }; // template DummyElectric

  template <typename T>
  const std::string DummyElectric<T>::tag = "DummyElectric";
  
  template <typename T> class DummyEx: public DummyElectric<T>
  {
  }; // template DummyEx

  template <typename T> class DummyEy: public DummyElectric<T>
  {
  }; // template DummyEy

  template <typename T> class DummyEz: public DummyElectric<T>
  {
  }; // template DummyEz

  template <typename T> 
  class DummyMagnetic: public MaterialMagnetic<T>
  {
  public:
    const std::string& 
    name() const
    {
      return DummyMagnetic<T>::tag;
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

      const auto& dummy_param = *static_cast<const DummyMagneticParam<T>*>(pm_param_ptr);
      
      idx_list.push_back(index);
      param_list.push_back(dummy_param);
      
      return this;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm_ptr)
    {
      auto dummy_ptr = static_cast<const DummyMagnetic<T>*>(pm_ptr);
      std::copy(dummy_ptr->idx_list.begin(), dummy_ptr->idx_list.end(), std::back_inserter(idx_list));
      std::copy(dummy_ptr->param_list.begin(), dummy_ptr->param_list.end(), std::back_inserter(param_list));
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
    using MaterialMagnetic<T>::position;
    using MaterialMagnetic<T>::idx_list;
    std::vector<DummyMagneticParam<T> > param_list;

  private:
    static const std::string tag; // "DummyMagnetic"
  }; // template DummyMagnetic

  template <typename T>
  const std::string DummyMagnetic<T>::tag = "DummyMagnetic";

  template <typename T> 
  class DummyHx: public DummyMagnetic<T>
  {
  }; // template DummyHx

  template <typename T> 
  class DummyHy: public DummyMagnetic<T>
  {
  }; // template DummyHy

  template <typename T> 
  class DummyHz: public DummyMagnetic<T>
  {
  }; // template DummyHz
} // namespace gmes

#endif // PW_DUMMY_HH_
