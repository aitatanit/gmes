#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

// #include <algorithm>
#include <array>
// #include <functional>
// #include <iterator>
// #include <utility>
// #include <vector>
#include <map>

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
  
  struct ltidx
  {
    bool
    operator()(const Index3& l, const Index3& r) const
    {
      return l[0] < r[0] || (l[0] == r[0] && (l[1] < r[1] || (l[1] == r[1] && l[2] < r[2])));
    }
  }; // struct ltidx

  typedef std::pair<Index3, PwMaterialParam*> IdxParamType;
  typedef std::map<Index3, PwMaterialParam*> IdxParamCont;
  // typedef std::vector<IdxParamType> IdxParamCont;

  // struct find_idx: std::unary_function<IdxParamType, bool> {
  //   Index3 idx;
  //   find_idx(const Index3& idx): idx(idx) {}

  //   bool operator() (const IdxParamType& ip) const {
  //     return ip.first == idx;
  //   }
  // }; // struct find_idx
  
  template<class T> 
  class PwMaterial 
  {
  public:
    virtual
    ~PwMaterial() {}

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
    
    IdxParamCont::const_iterator
    find(const Index3& idx) const
    {
      // auto it = std::find_if(param.begin(), param.end(), find_idx(idx));
      auto it = param.find(idx);
      return it;
    }

    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm)
    {
      // std::copy(pm->param.begin(), pm->param.end(), std::back_inserter(param));
      for (auto v: pm->param) {
        attach(v.first.data(), v.first.size(), v.second);
      }
      return this;
    }

    IdxParamCont::size_type
    idx_size() const
    {
      return param.size();
    }

    // IdxParamCont::const_iterator
    // begin() const
    // {
    //   return param.begin();
    // }
    
    // IdxParamCont::const_iterator
    // end() const
    // {
    //   return param.end();
    // }

  protected:
    IdxParamCont param;
  }; // template PwMaterial

  template <typename T> class MaterialElectric: public PwMaterial<T> 
  {
  public:
    double 
    get_eps_inf(const int* const idx, int idx_size) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());
      
      auto it = find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<ElectricParam<T>*>(it->second)->eps_inf;
    }

    using PwMaterial<T>::find;

  protected:
    using PwMaterial<T>::param;
  }; // template MaterialElectric

  template <typename T> class MaterialMagnetic: public PwMaterial<T> 
  {
  public:
    double 
    get_mu_inf(const int* const idx, int idx_size) const
    {
      Index3 index;
      std::copy(idx, idx + idx_size, index.begin());

      auto it = find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<MagneticParam<T>*>(it->second)->mu_inf;
    }

    using PwMaterial<T>::find;

  protected:
    using PwMaterial<T>::param;
  }; // template MaterialMagnetic
} // namespace gmes

#endif // PW_MATERIAL_HH_
