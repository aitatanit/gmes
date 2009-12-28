#ifndef PW_MATERIAL_HH_
#define PW_MATERIAL_HH_

#include <stdexcept>
#include "constants.hh"

namespace gmes {
template<class T> class PwMaterial {
public:
	// constructors & destructor
	PwMaterial(const int * const idx, int size)
	{
		if (size != 3)
		{
			throw std::domain_error("length of the array index must be 3.");
		}
		else
		{
			i = idx[0];
			j = idx[1];
			k = idx[2];
		}
	}

	virtual ~PwMaterial() {
	}

	int get_i() const {
		return i;
	}
	int get_j() const {
		return j;
	}
	int get_k() const {
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

	virtual void
	update(T * const inplace_field, int inplace_dim1, int inplace_dim2, int inplace_dim3,
			const T * const in_field1, int in1_dim1, int in1_dim2, int in1_dim3,
			const T * const in_field2, int in2_dim1, int in2_dim2, int in2_dim3,
			double d1, double d2, double dt, double t) = 0;

protected:
	int i, j, k;
};

template <class T> class MaterialElectric: public PwMaterial<T> {
public:
	MaterialElectric(const int * const idx, int size) :
		PwMaterial<T>(idx, size) {
	}
	virtual ~MaterialElectric() {
	}

	virtual double get_epsilon() = 0;
	virtual void set_epsilon(double epsilon_r) = 0;
};

template <class T> class MaterialMagnetic: public PwMaterial<T> {
public:
	MaterialMagnetic(const int * const idx, int size) :
		PwMaterial<T>(idx, size) {
	}
	virtual ~MaterialMagnetic() {
	}

	virtual double get_mu() = 0;
	virtual void set_mu(double mu_r) = 0;
};

}

#endif /*PW_MATERIAL_HH_*/
