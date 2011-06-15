#ifndef PW_DUMMY_HH_
#define PW_DUMMY_HH_

#include "pw_material.hh"

namespace gmes
{
  template <typename T> class DummyElectric: public MaterialElectric<T>
  {
  public:
    DummyElectric(double epsilon):
      eps(epsilon)
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

    void update(T * const inplace_field, 
		int inplace_dim1, int inplace_dim2, int inplace_dim3,
		const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
		const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
		double d1, double d2, double dt, double n, int i, int j, int k)
    {
      return;
    }

  protected:
    double eps;
  };

  template <typename T> class DummyEx: public DummyElectric<T>
  {
  public:
    DummyEx(double epsilon = 1):
      DummyElectric<T>(epsilon)
    {
    }
  };

  template <typename T> class DummyEy: public DummyElectric<T>
  {
  public:
    DummyEy(double epsilon = 1):
      DummyElectric<T>(epsilon)
    {
    }
  };

  template <typename T> class DummyEz: public DummyElectric<T>
  {
  public:
    DummyEz(double epsilon = 1):
      DummyElectric<T>(epsilon)
    {
    }
  };

  template <typename T> class DummyMagnetic: public MaterialMagnetic<T>
  {
  public:
    DummyMagnetic(double mu):
      mu(mu)
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

    void update(T * const inplace_field, 
		int inplace_dim1, int inplace_dim2, int inplace_dim3,
		const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
		const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
		double d1, double d2, double dt, double n, int i, int j, int k)
    {
      return;
    }

  protected:
    double mu;
  };

  template <typename T> class DummyHx: public DummyMagnetic<T>
  {
  public:
    DummyHx(double mu = 1):
      DummyMagnetic<T>(mu)
    {
    }
  };

  template <typename T> class DummyHy: public DummyMagnetic<T>
  {
  public:
    DummyHy(double mu = 1):
      DummyMagnetic<T>(mu)
    {
    }
  };

  template <typename T> class DummyHz: public DummyMagnetic<T>
  {
  public:
    DummyHz(double mu = 1):
      DummyMagnetic<T>(mu)
    {
    }
  };
}

#endif /*PW_DUMMY_HH_*/
