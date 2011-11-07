#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

#include <array>
#include <numeric>
#include <map>

namespace gmes 
{
  struct ltidx
  {
    bool 
    operator()(const std::array<int, 3>& l, const std::array<int, 3>& r) const
    {
      return (l[0] < r[0]) || (l[0] == r[0] && l[1] < r[1]) || (l[0] == r[0] && l[1] == r[1] && l[2] < r[2]);
    }
  };
  
  struct PwMaterialParam
  {
  };
  
  template <typename T> 
  struct ElectricParam: public PwMaterialParam
  {
    double eps_inf;
  };
  
  template <typename T> 
  struct MagneticParam: public PwMaterialParam
  {
    double mu_inf;
  };

  typedef std::map<std::array<int, 3>, PwMaterialParam*, ltidx> MapType;

  template<class T> 
  class PwMaterial 
  {
  public:
    virtual
    ~PwMaterial()
    {
    }

    virtual PwMaterial<T>*
    attach(const int* const idx, int idx_size,
	   const PwMaterialParam* const parameter) = 0;
    
    void
    update_all(T* const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T* const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T* const in_field2, 
	       int in2_dim1, int in2_dim2, int in2_dim3,
	       double d1, double d2, double dt, double n)
    {
      for(MapType::const_iterator iter = param.begin(); 
	  iter != param.end(); iter++) {
	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
	       in_field1, in1_dim1, in1_dim2, in1_dim3,
	       in_field2, in2_dim1, in2_dim2, in2_dim3,
	       d1, d2, dt, n, 
	       iter->first.data(), iter->first.size(), iter->second);
      }
    }
    
    bool
    has_it(const int* const idx, int idx_size) const
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      return param.find(index) != param.end();
    }
    
    PwMaterial<T>*
    merge(const PwMaterial<T>* const pm)
    {
      for (MapType::const_iterator it = pm->param.begin();
      	   it != pm->param.end(); it++) {
      	attach(it->first.data(), it->first.size(), it->second);
      }
      return this;
    }

    MapType::size_type
    idx_size() const
    {
      return param.size();
    }

    MapType::const_iterator
    begin() const
    {
      return param.begin();
    }
    
    MapType::const_iterator
    end() const
    {
      return param.end();
    }

    // typedef typename MapType::const_iterator const_iterator;

  protected:
    virtual void
    update(T* const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T* const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T* const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int* const idx, int idx_size,
	   PwMaterialParam* const parameter) = 0;

    MapType param;
  };

  template <typename T> class MaterialElectric: public PwMaterial<T> 
  {
  public:
    ~MaterialElectric()
    {
    }

    double 
    get_eps_inf(const int* const idx, int idx_size) const
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<ElectricParam<T>*>(it->second)->eps_inf;
    }

  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class MaterialMagnetic: public PwMaterial<T> 
  {
  public:
    ~MaterialMagnetic() 
    {
    }

    double 
    get_mu_inf(const int* const idx, int idx_size) const
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<MagneticParam<T>*>(it->second)->mu_inf;
    }

  protected:
    using PwMaterial<T>::param;
  };
}

#endif /*PW_MATERIAL_HH_*/
