#ifndef PW_DUMMY_REAL_HH_
#define PW_DUMMY_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class DummyElectricReal: public MaterialElectricReal
{
public:
	DummyElectricReal(const int * const idx, int size, double epsilon_r) :
		MaterialElectricReal(idx, size), epsilon(epsilon_r * epsilon0)
	{
	}

	double get_epsilon()
	{
		return epsilon;
	}

	void set_epsilon(double epsilon_r)
	{
		epsilon = epsilon_r * epsilon0;
	}

	void update(double * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const double * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const double * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double t);

protected:
	double epsilon;
};

class DummyExReal: public DummyElectricReal
{
public:
	DummyExReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}
};

class DummyEyReal: public DummyElectricReal
{
public:
	DummyEyReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}
};

class DummyEzReal: public DummyElectricReal
{
public:
	DummyEzReal(const int * const idx, int size, double epsilon_r = 1) :
		DummyElectricReal(idx, size, epsilon_r)
	{
	}
};

class DummyMagneticReal: public MaterialMagneticReal
{
public:
	DummyMagneticReal(const int * const idx, int size, double mu_r) :
		MaterialMagneticReal(idx, size), mu(mu_r * mu0)
	{
	}

	double get_mu()
	{
		return mu;
	}

	void set_mu(double mu_r)
	{
		mu = mu_r * mu0;
	}

	void update(double * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
				const double * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
				const double * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
				double d1, double d2, double dt, double t);

protected:
	double mu;
};

class DummyHxReal: public DummyMagneticReal
{
public:
	DummyHxReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}
};

class DummyHyReal: public DummyMagneticReal
{
public:
	DummyHyReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}
};

class DummyHzReal: public DummyMagneticReal
{
public:
	DummyHzReal(const int * const idx, int size, double mu_r = 1) :
		DummyMagneticReal(idx, size, mu_r)
	{
	}
};
}

#endif /*PW_DUMMY_REAL_HH_*/
