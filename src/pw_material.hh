#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

#include <algorithm>
#include <array>
#include <iterator>
#include <functional>
#include <utility>
#include <vector>

namespace gmes 
{
  struct PwMaterialParam
  {
  }; // struct PwMaterialParam
  
  template <typename T> 
  struct ElectricParam: public PwMaterialParam
  {
    double eps_inf;
  }; // template ElectricParam
  
  template <typename T> 
  struct MagneticParam: public PwMaterialParam
  {
    double mu_inf;
  }; // template MagneticParam

  typedef std::array<int, 3> Index3;
  typedef std::vector<Index3> IdxCnt;

  template <typename T> 
  class PwMaterial 
  {
  public:
    virtual
    ~PwMaterial() {}

    virtual const std::string& name() const = 0;

    // TODO: just copy PwMaterialParam*.
    virtual PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const parameter) = 0;
    
    virtual void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n) = 0;
    
    IdxCnt::const_iterator
    find(const Index3& idx) const
    {
      auto it = std::find(idx_list.begin(), idx_list.end(), idx);
      return it;
    }

    virtual PwMaterial<T>*
    merge(const PwMaterial<T>* const pm) = 0;

    IdxCnt::size_type
    idx_size() const
    {
      return idx_list.size();
    }

  protected:
    int
    position(const Index3& idx) const
    {
      auto it = find(idx);
      
      size_t pos = 0;
      if (it == idx_list.end())
	return pos - 1;
      else {
	pos = std::distance(idx_list.begin(), it);
	return pos;
      }
    }
    
    IdxCnt idx_list;
  }; // template PwMaterial

  template <typename T> 
  class MaterialElectric: public PwMaterial<T> 
  {
  public:
    virtual double 
    get_eps_inf(const int* const idx, int idx_size) const = 0;

    using PwMaterial<T>::find;

  protected:
    using PwMaterial<T>::position;
    using PwMaterial<T>::idx_list;
  }; // template MaterialElectric

  template <typename T> 
  class MaterialMagnetic: public PwMaterial<T> 
  {
  public:
    virtual double 
    get_mu_inf(const int* const idx, int idx_size) const = 0;

    using PwMaterial<T>::find;

  protected:
    using PwMaterial<T>::position;
    using PwMaterial<T>::idx_list;
  }; // template MaterialMagnetic
} // namespace gmes

#endif // PW_MATERIAL_HH_
