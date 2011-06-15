#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

namespace gmes {
  template<class T> class PwMaterial {
  public:
    virtual ~PwMaterial() {
    }

    virtual void
    update(T * const inplace_field, 
	   int inplace_dim1, int inplace_dim2, int inplace_dim3,
	   const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
	   const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
	   double d1, double d2, double dt, double n, int i, int j, int k) = 0;
  };

  template <typename T> class MaterialElectric: public PwMaterial<T> {
  public:
    virtual ~MaterialElectric() {
    }

    virtual double get_epsilon() const = 0;
    virtual void set_epsilon(double epsilon) = 0;
  };

  template <typename T> class MaterialMagnetic: public PwMaterial<T> {
  public:
    virtual ~MaterialMagnetic() {
    }

    virtual double get_mu() const = 0;
    virtual void set_mu(double mu) = 0;
  };

}

#endif /*PW_MATERIAL_HH_*/
