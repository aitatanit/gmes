#ifndef PW_ZERO_REAL_HH_
#define PW_ZERO_REAL_HH_

#include "pw_material_real.hh"
#include "constants.hh"

namespace gmes
{
class ZeroElectricReal: public MaterialElectricReal
{
public:
	ZeroElectricReal(const int * const idx, int size, double epsilon_r) :
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

class ZeroExReal: public ZeroElectricReal
{
public:
	ZeroExReal(const int * const idx, int size, double epsilon_r = 1) :
			ZeroElectricReal(idx, size, epsilon_r)
	{
	}
};

class ZeroEyReal: public ZeroElectricReal
{
public:
	ZeroEyReal(const int * const idx, int size, double epsilon_r = 1) :
			ZeroElectricReal(idx, size, epsilon_r)
	{
	}
};

class ZeroEzReal: public ZeroElectricReal
{
public:
	ZeroEzReal(const int * const idx, int size, double epsilon_r = 1) :
			ZeroElectricReal(idx, size, epsilon_r)
	{
	}
};

class ZeroMagneticReal: public MaterialMagneticReal
{
public:
	ZeroMagneticReal(const int * const idx, int size, double mu_r) :
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

class ZeroHxReal: public ZeroMagneticReal
{
public:
	ZeroHxReal(const int * const idx, int size, double mu_r = 1) :
			ZeroMagneticReal(idx, size, mu_r)
	{
	}
};

class ZeroHyReal: public ZeroMagneticReal
{
public:
	ZeroHyReal(const int * const idx, int size, double mu_r = 1) :
			ZeroMagneticReal(idx, size, mu_r)
	{
	}
};

class ZeroHzReal: public ZeroMagneticReal
{
public:
	ZeroHzReal(const int * const idx, int size, double mu_r = 1) :
			ZeroMagneticReal(idx, size, mu_r)
	{
	}
};
}

#endif /*PW_ZERO_REAL_HH_*/
