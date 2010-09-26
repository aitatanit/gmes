#ifndef PW_DUMMY_HH_
#define PW_DUMMY_HH_

#include "pw_material.hh"

namespace gmes
{
template <typename T> class DummyElectric: public MaterialElectric<T>
{
public:
	DummyElectric(const int * const idx, int size, double epsilon) :
		MaterialElectric<T>(idx, size), epsilon(epsilon)
	{
	}

	double get_epsilon() const
	{
		return epsilon;
	}

	void set_epsilon(double epsilon)
	{
		epsilon = epsilon;
	}

	void update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double n)
	{
		return;
	}

protected:
	double epsilon;
};

template <typename T> class DummyEx: public DummyElectric<T>
{
public:
	DummyEx(const int * const idx, int size, double epsilon = 1) :
		DummyElectric<T>(idx, size, epsilon)
	{
	}
};

template <typename T> class DummyEy: public DummyElectric<T>
{
public:
	DummyEy(const int * const idx, int size, double epsilon = 1) :
		DummyElectric<T>(idx, size, epsilon)
	{
	}
};

template <typename T> class DummyEz: public DummyElectric<T>
{
public:
	DummyEz(const int * const idx, int size, double epsilon = 1) :
		DummyElectric<T>(idx, size, epsilon)
	{
	}
};

template <typename T> class DummyMagnetic: public MaterialMagnetic<T>
{
public:
	DummyMagnetic(const int * const idx, int size, double mu) :
		MaterialMagnetic<T>(idx, size), mu(mu)
	{
	}

	double get_mu() const
	{
		return mu;
	}

	void set_mu(double mu)
	{
		mu = mu;
	}

	void update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double n)
	{
		return;
	}

protected:
	double mu;
};

template <typename T> class DummyHx: public DummyMagnetic<T>
{
public:
	DummyHx(const int * const idx, int size, double mu = 1) :
		DummyMagnetic<T>(idx, size, mu)
	{
	}
};

template <typename T> class DummyHy: public DummyMagnetic<T>
{
public:
	DummyHy(const int * const idx, int size, double mu = 1) :
		DummyMagnetic<T>(idx, size, mu)
	{
	}
};

template <typename T> class DummyHz: public DummyMagnetic<T>
{
public:
	DummyHz(const int * const idx, int size, double mu = 1) :
		DummyMagnetic<T>(idx, size, mu)
	{
	}
};
}

#endif /*PW_DUMMY_HH_*/
