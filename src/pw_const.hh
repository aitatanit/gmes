#ifndef PW_CONST_HH_
#define PW_CONST_HH_

#include "pw_material.hh"

#define inplace_field(i,j,k) inplace_field[((this->i)*inplace_dim2+(this->j))*inplace_dim3+(this->k)]

namespace gmes
{
template <typename T> class ConstElectric: public MaterialElectric<T>
{
public:
	ConstElectric(const int * const idx, int size, double epsilon, T value) :
		MaterialElectric<T>(idx, size), epsilon(epsilon), value(value)
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

	T get_value() const
	{
		return value;
	}

	void set_value(T value)
	{
		this->value = value;
	}

	void update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double t)
	{
		inplace_field(i,j,k) = value;
	}

protected:
	double epsilon;
	T value;
};

template <typename T> class ConstEx: public ConstElectric<T>
{
public:
	ConstEx(const int * const idx, int size, double epsilon = 1, const T& value = 0.) :
			ConstElectric<T>(idx, size, epsilon, value)
	{
	}
};

template <typename T> class ConstEy: public ConstElectric<T>
{
public:
	ConstEy(const int * const idx, int size, double epsilon = 1, const T& value = 0.) :
			ConstElectric<T>(idx, size, epsilon, value)
	{
	}
};

template <typename T> class ConstEz: public ConstElectric<T>
{
public:
	ConstEz(const int * const idx, int size, double epsilon = 1, const T& value = 0.) :
			ConstElectric<T>(idx, size, epsilon, value)
	{
	}
};

template <typename T> class ConstMagnetic: public MaterialMagnetic<T>
{
public:
	ConstMagnetic(const int * const idx, int size, double mu, T value) :
		MaterialMagnetic<T>(idx, size), mu(mu), value(value)
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

	T get_value() const
	{
		return value;
	}

	void set_value(T value)
	{
		this->value = value;
	}

	void update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double t)

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
	ConstHx(const int * const idx, int size, double mu = 1, const T& value = 0.) :
			ConstMagnetic<T>(idx, size, mu, value)
	{
	}
};

template <typename T> class ConstHy: public ConstMagnetic<T>
{
public:
	ConstHy(const int * const idx, int size, double mu = 1, const T& value = 0.) :
			ConstMagnetic<T>(idx, size, mu, value)
	{
	}
};

template <typename T> class ConstHz: public ConstMagnetic<T>
{
public:
	ConstHz(const int * const idx, int size, double mu = 1, const T& value = 0.) :
			ConstMagnetic<T>(idx, size, mu, value)
	{
	}
};
}

#undef inplace_field

#endif /*PW_CONST_HH_*/
