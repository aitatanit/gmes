#ifndef PW_MATERIAL_REAL_HH_
#define PW_MATERIAL_REAL_HH_

#include "constants.hh"

namespace gmes {
class PwMaterialReal {
public:
	// constructors & destructor
	PwMaterialReal(const int * const idx, int size);
	virtual ~PwMaterialReal() {
	}

	int get_i() {
		return i;
	}
	int get_j() {
		return j;
	}
	int get_k() {
		return k;
	}
	void set_i(int i) {
		this->i = i;
	}
	void set_j(int j) {
		this->j = j;
	}
	void set_k(int k) {
		this->k = k;
	}

	virtual void update(double * const inplace_field, int inplace_dim1,
			int inplace_dim2, int inplace_dim3, const double * const in_field1,
			int in1_dim1, int in1_dim2, int in1_dim3,
			const double * const in_field2, int in2_dim1, int in2_dim2,
			int in2_dim3, double, double, double, double) = 0;

protected:
	int i, j, k;
};

class MaterialElectricReal: public PwMaterialReal {
public:
	MaterialElectricReal(const int * const idx, int size) :
		PwMaterialReal(idx, size) {
	}

	virtual ~MaterialElectricReal() {
	}

	virtual double get_epsilon() = 0;
	virtual void set_epsilon(double epsilon_r) = 0;
};

class MaterialMagneticReal: public PwMaterialReal {
public:
	MaterialMagneticReal(const int * const idx, int size) :
		PwMaterialReal(idx, size) {
	}

	virtual ~MaterialMagneticReal() {
	}

	virtual double get_mu() = 0;
	virtual void set_mu(double mu_r) = 0;
};

}

#endif /*PW_MATERIAL_REAL_HH_*/
