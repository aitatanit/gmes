#ifndef PW_CONST_HH_
#define PW_CONST_HH_

#include "pw_material.hh"

#define inplace_field(i,j,k) inplace_field[((i)*inplace_dim2+(j))*inplace_dim3+(k)]

namespace gmes
{
  template <typename T> class ConstElectric: public MaterialElectric<T>
  {
  public:
    ConstElectric(double epsilon, T value):
      eps(epsilon), value(value)
    {
    }

    double get_epsilon() const
    {
      return eps;
    }

    void set_epsilon(double epsilon)
    {
      eps = epsilon;
    }

    T get_value() const
    {
      return value;
    }

    void set_value(T value)
    {
      this->value = value;
    }

    void update(T * const inplace_field, 
		int inplace_dim1, int inplace_dim2, int inplace_dim3,
		const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
		const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
		double d1, double d2, double dt, double n, int i, int j, int k)
    {
      inplace_field(i,j,k) = value;
    }

  protected:
    double eps;
    T value;
  };

  template <typename T> class ConstEx: public ConstElectric<T>
  {
  public:
    ConstEx(double epsilon = 1, const T& value = 0):
      ConstElectric<T>(epsilon, value)
    {
    }
  };

  template <typename T> class ConstEy: public ConstElectric<T>
  {
  public:
    ConstEy(double epsilon = 1, const T& value = 0):
      ConstElectric<T>(epsilon, value)
    {
    }
  };

  template <typename T> class ConstEz: public ConstElectric<T>
  {
  public:
    ConstEz(double epsilon = 1, const T& value = 0):
      ConstElectric<T>(epsilon, value)
    {
    }
  };

  template <typename T> class ConstMagnetic: public MaterialMagnetic<T>
  {
  public:
    ConstMagnetic(double mu, T value):
      mu(mu), value(value)
    {
    }

    double get_mu() const
    {
      return mu;
    }

    void set_mu(double mu)
    {
      this->mu = mu;
    }

    T get_value() const
    {
      return value;
    }

    void set_value(T value)
    {
      this->value = value;
    }

    void update(T * const inplace_field, 
		int inplace_dim1, int inplace_dim2, int inplace_dim3,
		const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
		const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
		double d1, double d2, double dt, double n, int i, int j, int k)

    {
      inplace_field(i,j,k) = value;
    }

  protected:
    double mu;
    T value;
  };

  template <typename T> class ConstHx: public ConstMagnetic<T>
  {
  public:
    ConstHx(double mu = 1, const T& value = 0):
      ConstMagnetic<T>(mu, value)
    {
    }
  };

  template <typename T> class ConstHy: public ConstMagnetic<T>
  {
  public:
    ConstHy(double mu = 1, const T& value = 0):
      ConstMagnetic<T>(mu, value)
    {
    }
  };

  template <typename T> class ConstHz: public ConstMagnetic<T>
  {
  public:
    ConstHz(double mu = 1, const T& value = 0):
      ConstMagnetic<T>(mu, value)
    {
    }
  };
}

#undef inplace_field

#endif /*PW_CONST_HH_*/
