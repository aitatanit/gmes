#ifndef POINTWISE_CPML_HH_
#define POINTWISE_CPML_HH_

#include "pointwise_material.hh"

namespace gmes 
{
    class CPMLElectric: public PointwiseMaterial
    {
    public:
	CPMLElectric(const int * const idx, int size, double epsilon_r, 
		     double b1, double b2, double c1, double c2, double kappa1, double kappa2);
	double get_epsilon() { return epsilon; }
	void set_epsilon(double epsilon_r);

    protected:
	double epsilon;
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	double psi1, psi2;
    };


    class CPMLEx: public CPMLElectric
    {
    public:
	CPMLEx(const int * const idx, int size, double epsilon_r, 
	       double by, double bz, double cy, double cz, double kappay, double kappaz)
	    : CPMLElectric(idx, size, epsilon_r, by, bz, cy, cz, kappay, kappaz) {}

	void update(double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		    const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		    const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		    double dt, double dy, double dz);
    };


    class CPMLEy: public CPMLElectric
    {
    public:
	CPMLEy(const int * const idx, int size,  double epsilon_r, 
	       double bz, double bx, double cz, double cx, double kappaz, double kappax)
	    : CPMLElectric(idx, size, epsilon_r, bz, bx, cz, cx, kappaz, kappax) {}

	void update(double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		    const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		    const double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		    double dt, double dz, double dx);
    };


    class CPMLEz: public CPMLElectric
    {
    public:
	CPMLEz(const int * const idx, int size, double epsilon_r, 
	       double bx, double by, double cx, double cy, double kappax, double kappay)
	    : CPMLElectric(idx, size, epsilon_r, bx, by, cx, cy, kappax, kappay) {}

	void update(double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		    const double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		    const double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		    double dt, double dx, double dy);
    };


    class CPMLMagnetic: public PointwiseMaterial
    {
    public:
	CPMLMagnetic(const int * const idx, int size, double mu_r, 
		     double b1, double b2, double c1, double c2, double kappa1, double kappa2);
	double get_mu() { return mu; }
	void set_mu(double mu_r);

    protected:
	double mu;
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	double psi1, psi2;
    };


    class CPMLHx: public CPMLMagnetic
    {
    public:
	CPMLHx(const int * const idx, int size, double mu_r, 
	       double by, double bz, double cy, double cz, double kappay, double kappaz)
	    : CPMLMagnetic(idx, size, mu_r, by, bz, cy, cz, kappay, kappaz) {}

	void update(double * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
		    const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		    const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,	    
		    double dt, double dy, double dz);
    };


    class CPMLHy: public CPMLMagnetic
    {
    public:
	CPMLHy(const int * const idx, int size, double mu_r, 
	       double bz, double bx, double cz, double cx, double kappaz, double kappax)
	    : CPMLMagnetic(idx, size, mu_r, bz, bx, cz, cx, kappaz, kappax) {}

	void update(double * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
		    const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		    const double * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
		    double dt, double dz, double dx);
    };


    class CPMLHz: public CPMLMagnetic
    {
    public:
	CPMLHz(const int * const idx, int size, double mu_r, 
	       double bx, double by, double cx, double cy, double kappax, double kappay)
	    : CPMLMagnetic(idx, size, mu_r, bx, by, cx, cy, kappax, kappay) {}

	void update(double * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
		    const double * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
		    const double * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
		    double dt, double dx, double dy);
    };
}

#endif /*POINTWISE_CPML_HH_*/
