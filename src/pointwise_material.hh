#ifndef POINTWISE_MATERIAL_HH_
#define POINTWISE_MATERIAL_HH_

namespace gmes
{
    class PointwiseMaterial
    {
    public:
	// constructors & destructors
	PointwiseMaterial(const int * const idx, int size);
	virtual ~PointwiseMaterial() {};
	
	int get_i() { return i; }
	int get_j() { return j; }
	int get_k() { return k; }
	void set_i(int i) { this->i = i; }
	void set_j(int j) { this->j = j; }
	void set_k(int k) { this->k = k; }

	virtual void update(double * const inplace_field, 
			    int inplace_dim1, int inplace_dim2, int inplace_dim3,
			    const double * const in_field1, 
			    int in1_dim1, int in1_dim2, int in1_dim3,
			    const double * const in_field2, 
			    int in2_dim1, int in2_dim2, int in2_dim3,
			    double, double, double) = 0;

    protected:
	int i, j, k;
    };
}

#endif /*POINTWISE_MATERIAL_HH_*/
