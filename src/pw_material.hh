#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

#include <unordered_map>

namespace gmes 
{
  struct eqidx{
    bool operator()(const int idx1[3], const int idx2[3]) const {
      return idx1[0] == idx2[0] && idx1[1] == idx2[1] && idx1[2] == idx2[2];
    }
  };
  
  struct PwMaterialParam 
  {
  };

  struct ElectricParam: public PwMaterialParam 
  {
    double eps;
  };

  struct MagneticParam: public PwMaterialParam 
  {
    double mu;
  };

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
    
    virtual void
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, 
	   const int idx[3], int idx_size,
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
      for(std::unordered_map<const int[3], PwMaterialParam *, eqidx>::const_iterator iter = param.begin(); iter != param.end(); iter++) {
	update(inplace_field, inplace_dim1, inplace_dim2, inplace_dim3,
	       in_field1, in1_dim1, in1_dim2, in1_dim3,
	       in_field2, in2_dim1, in2_dim2, in2_dim3,
	       d1, d2, dt, n, iter->first, iter->second);
      }
    }
    
    bool
    has_it(const int idx[3], int idx_size)
    {
      return param.find(const_cast<int[3]>(idx)) != param.end();
    }
    
  protected:
    std::unordered_map<const int[3], PwMaterialParam *, eqidx> param;
  };

  template <typename T> class MaterialElectric: public PwMaterial<T> 
  {
  public:
    virtual 
    ~MaterialElectric()
    {
    }

    double 
    get_epsilon(const int idx[3], int idx_size) const
    {
      return static_cast<ElectricParam *>(param[idx])->eps;
    }

    double 
    set_epsilon(const int idx[3], int idx_size, double epsilon)
    {
      return static_cast<ElectricParam *>(param[idx]) = epsilon;
    }

  protected:
    using PwMaterial<T>::param;
  };

  template <typename T> class MaterialMagnetic: public PwMaterial<T> 
  {
  public:
    virtual 
    ~MaterialMagnetic() 
    {
    }

    double get_mu(const int idx[3], int idx_size) const
    {
      return static_cast<MagneticParam *>(param[idx])->mu;
    }

    double set_mu(const int idx[3], int idx_size, double mu)
    {
      return static_cast<MagneticParam *>(param[idx]) = mu;
    }

  protected:
    using PwMaterial<T>::param;
  };

  typedef std::unordered_map<const int[3], PwMaterialParam *, eqidx> MapType;
}

#endif /*PW_MATERIAL_HH_*/
