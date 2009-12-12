#ifndef PW_ONE_REAL_HH_
#define PW_ONE_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class OneElectricReal: public MaterialElectricReal
{
public:
	OneElectricReal(const int * const idx, int size, double epsilon_r) :
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

class OneExReal: public OneElectricReal
{
public:
	OneExReal(const int * const idx, int size, double epsilon_r = 1) :
			OneElectricReal(idx, size, epsilon_r)
	{
	}
};

class OneEyReal: public OneElectricReal
{
public:
	OneEyReal(const int * const idx, int size, double epsilon_r = 1) :
			OneElectricReal(idx, size, epsilon_r)
	{
	}
};

class OneEzReal: public OneElectricReal
{
public:
	OneEzReal(const int * const idx, int size, double epsilon_r = 1) :
			OneElectricReal(idx, size, epsilon_r)
	{
	}
};

class OneMagneticReal: public MaterialMagneticReal
{
public:
	OneMagneticReal(const int * const idx, int size, double mu_r) :
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

class OneHxReal: public OneMagneticReal
{
public:
	OneHxReal(const int * const idx, int size, double mu_r = 1) :
			OneMagneticReal(idx, size, mu_r)
	{
	}
};

class OneHyReal: public OneMagneticReal
{
public:
	OneHyReal(const int * const idx, int size, double mu_r = 1) :
			OneMagneticReal(idx, size, mu_r)
	{
	}
};

class OneHzReal: public OneMagneticReal
{
public:
	OneHzReal(const int * const idx, int size, double mu_r = 1) :
			OneMagneticReal(idx, size, mu_r)
	{
	}
};
}

#endif /*PW_ONE_REAL_HH_*/
