#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

#include <array>
#include <numeric>
#include <unordered_map>

namespace std
{
  template <>
  struct hash<array<int, 3> >: public unary_function<array<int, 3>, size_t>
  {
    size_t operator()(const array<int, 3>& idx) const
    {
      return std::accumulate(idx.begin(), idx.end(), 0);
    }
  };
}

namespace gmes 
{
  struct PwMaterialParam
  {
  };
  
  template <typename T> struct ElectricParam: public PwMaterialParam
  {
    double eps;
  };
  
  template <typename T> struct MagneticParam: public PwMaterialParam
  {
    double mu;
  };

  typedef std::unordered_map<std::array<int, 3>, PwMaterialParam *> MapType;

  template<class T> class PwMaterial 
  {
  public:
    virtual
    ~PwMaterial()
    {
    }

    virtual void 
    attach(const int idx[3], int idx_size,
	   const PwMaterialParam * const parameter) = 0;
    
    void
    update_all(T * const inplace_field,
	       int inplace_dim1, int inplace_dim2, int inplace_dim3,
	       const T * const in_field1, 
	       int in1_dim1, int in1_dim2, int in1_dim3,
	       const T * const in_field2, 
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
    has_it(const int idx[3], int idx_size)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      return param.find(index) != param.end();
    }
    
  protected:
    virtual void
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int idx[3], int idx_size,
	   PwMaterialParam * const parameter) = 0;

    MapType param;
  };

  template <typename T> class MaterialElectric: public PwMaterial<T> 
  {
  public:
    ~MaterialElectric()
    {
    }

    double 
    get_epsilon(const int idx[3], int idx_size) const
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<ElectricParam<T> *>(it->second)->eps;
    }

    double 
    set_epsilon(const int idx[3], int idx_size, double epsilon)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());
      
      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<ElectricParam<T> *>(it->second)->eps = epsilon;
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
    get_mu(const int idx[3], int idx_size) const
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<MagneticParam<T> *>(it->second)->mu;
    }

    double 
    set_mu(const int idx[3], int idx_size, double mu)
    {
      std::array<int, 3> index;
      std::copy(idx, idx + idx_size, index.begin());

      MapType::const_iterator it = param.find(index);
      if (it == param.end())
	return 0;
      else
	return static_cast<MagneticParam<T> *>(it->second)->mu = mu;
    }

  protected:
    using PwMaterial<T>::param;
  };
}

#endif /*PW_MATERIAL_HH_*/
