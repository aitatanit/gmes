/* This implementation is based on the following article.
 * S. Gedney, “Perfectly Matched Layer Absorbing Boundary Conditions,”
 * Computational Electrodynamics: The Finite-Difference Time-Domain Method,
 * Third Edition, A. Taflove and S.C. Hagness, eds., Artech House Publishers,
 *  2005, pp. 273-328.
 */

#ifndef PW_CPML_HH_
#define PW_CPML_HH_

#include "pw_material.hh"

#define ex(i,j,k) ex[((this->i)*ex_y_size+(this->j))*ex_z_size+(this->k)]
#define ey(i,j,k) ey[((this->i)*ey_y_size+(this->j))*ey_z_size+(this->k)]
#define ez(i,j,k) ez[((this->i)*ez_y_size+(this->j))*ez_z_size+(this->k)]
#define hx(i,j,k) hx[((this->i)*hx_y_size+(this->j))*hx_z_size+(this->k)]
#define hy(i,j,k) hy[((this->i)*hy_y_size+(this->j))*hy_z_size+(this->k)]
#define hz(i,j,k) hz[((this->i)*hz_y_size+(this->j))*hz_z_size+(this->k)]

namespace gmes
{
template <typename T> class CpmlElectric: public MaterialElectric<T>
{
public:
	CpmlElectric(const int * const idx, int size, double epsilon,
			double b1_in, double b2_in, double c1_in, double c2_in,
			double kappa1_in, double kappa2_in) :
		MaterialElectric<T>(idx, size), epsilon(epsilon),
		b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
		kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
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

protected:
	double epsilon;
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	T psi1, psi2;
};

template <typename T> class CpmlEx: public CpmlElectric<T>
{
public:
	CpmlEx(const int * const idx, int size, double epsilon, double by,
			double bz, double cy, double cz, double kappay, double kappaz) :
		CpmlElectric<T>(idx, size, epsilon, by, bz, cy, cz, kappay, kappaz)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (hz(i+1,j+1,k) - hz(i+1,j,k)) / dy;
		psi2 = b2 * psi2 + c2 * (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz;

		ex(i,j,k) += dt / epsilon * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy
				/ kappa1 - (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz / kappa2 + psi1
				- psi2);
	}

protected:
	using CpmlElectric<T>::epsilon;
	using CpmlElectric<T>::b1;
	using CpmlElectric<T>::b2;
	using CpmlElectric<T>::c1;
	using CpmlElectric<T>::c2;
	using CpmlElectric<T>::kappa1;
	using CpmlElectric<T>::kappa2;
	using CpmlElectric<T>::psi1;
	using CpmlElectric<T>::psi2;
};

template <typename T> class CpmlEy: public CpmlElectric<T>
{
public:
	CpmlEy(const int * const idx, int size, double epsilon, double bz,
			double bx, double cz, double cx, double kappaz, double kappax) :
		CpmlElectric<T>(idx, size, epsilon, bz, bx, cz, cx, kappaz, kappax)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (hx(i,j+1,k+1) - hx(i,j+1,k)) / dz;
		psi2 = b2 * psi2 + c2 * (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx;

		ey(i,j,k) += dt / epsilon * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz
				/ kappa1 - (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx / kappa2 + psi1
				- psi2);
	}

protected:
	using CpmlElectric<T>::epsilon;
	using CpmlElectric<T>::b1;
	using CpmlElectric<T>::b2;
	using CpmlElectric<T>::c1;
	using CpmlElectric<T>::c2;
	using CpmlElectric<T>::kappa1;
	using CpmlElectric<T>::kappa2;
	using CpmlElectric<T>::psi1;
	using CpmlElectric<T>::psi2;
};

template <typename T> class CpmlEz: public CpmlElectric<T>
{
public:
	CpmlEz(const int * const idx, int size, double epsilon, double bx,
			double by, double cx, double cy, double kappax, double kappay) :
		CpmlElectric<T>(idx, size, epsilon, bx, by, cx, cy, kappax, kappay)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (hy(i+1,j,k+1) - hy(i,j,k+1)) / dx;
		psi2 = b2 * psi2 + c2 * (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy;

		ez(i,j,k) += dt / epsilon * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx
				/ kappa1 - (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy / kappa2 + psi1
				- psi2);
	}

protected:
	using CpmlElectric<T>::epsilon;
	using CpmlElectric<T>::b1;
	using CpmlElectric<T>::b2;
	using CpmlElectric<T>::c1;
	using CpmlElectric<T>::c2;
	using CpmlElectric<T>::kappa1;
	using CpmlElectric<T>::kappa2;
	using CpmlElectric<T>::psi1;
	using CpmlElectric<T>::psi2;
};

template <typename T> class CpmlMagnetic: public MaterialMagnetic<T>
{
public:
	CpmlMagnetic(const int * const idx, int size, double mu, double b1_in,
			double b2_in, double c1_in, double c2_in, double kappa1_in,
			double kappa2_in) :
		MaterialMagnetic<T>(idx, size), mu(mu),
		b1(b1_in), b2(b2_in), c1(c1_in), c2(c2_in),
		kappa1(kappa1_in), kappa2(kappa2_in), psi1(0), psi2(0)
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

protected:
	double mu;
	double b1, b2;
	double c1, c2;
	double kappa1, kappa2;
	T psi1, psi2;
};

template <typename T> class CpmlHx: public CpmlMagnetic<T>
{
public:
	CpmlHx(const int * const idx, int size, double mu, double by, double bz,
			double cy, double cz, double kappay, double kappaz) :
		CpmlMagnetic<T>(idx, size, mu, by, bz, cy, cz, kappay, kappaz)
	{
	}

	void update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dy, double dz, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy;
		psi2 = b2 * psi2 + c2 * (ey(i,j-1,k) - ey(i,j-1,k-1)) / dz;

		hx(i,j,k) -= dt / mu
				* ((ez(i,j,k-1) - ez(i,j-1,k-1)) / dy / kappa1 - (ey(i,j-1,k)
						- ey(i,j-1,k-1)) / dz / kappa2 + psi1 - psi2);
	}

protected:
	using CpmlMagnetic<T>::mu;
	using CpmlMagnetic<T>::b1;
	using CpmlMagnetic<T>::b2;
	using CpmlMagnetic<T>::c1;
	using CpmlMagnetic<T>::c2;
	using CpmlMagnetic<T>::kappa1;
	using CpmlMagnetic<T>::kappa2;
	using CpmlMagnetic<T>::psi1;
	using CpmlMagnetic<T>::psi2;
};

template <typename T> class CpmlHy: public CpmlMagnetic<T>
{
public:
	CpmlHy(const int * const idx, int size, double mu, double bz, double bx,
			double cz, double cx, double kappaz, double kappax) :
		CpmlMagnetic<T>(idx, size, mu, bz, bx, cz, cx, kappaz, kappax)
	{
	}

	void update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dz, double dx, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz;
		psi2 = b2 * psi2 + c2 * (ez(i,j,k-1) - ez(i-1,j,k-1)) / dx;

		hy(i,j,k) -= dt / mu
				* ((ex(i-1,j,k) - ex(i-1,j,k-1)) / dz / kappa1 - (ez(i,j,k-1)
						- ez(i-1,j,k-1)) / dx / kappa2 + psi1 - psi2);
	}

protected:
	using CpmlMagnetic<T>::mu;
	using CpmlMagnetic<T>::b1;
	using CpmlMagnetic<T>::b2;
	using CpmlMagnetic<T>::c1;
	using CpmlMagnetic<T>::c2;
	using CpmlMagnetic<T>::kappa1;
	using CpmlMagnetic<T>::kappa2;
	using CpmlMagnetic<T>::psi1;
	using CpmlMagnetic<T>::psi2;
};

template <typename T> class CpmlHz: public CpmlMagnetic<T>
{
public:
	CpmlHz(const int * const idx, int size, double mu, double bx, double by,
			double cx, double cy, double kappax, double kappay) :
		CpmlMagnetic<T>(idx, size, mu, bx, by, cx, cy, kappax, kappay)
	{
	}

	void update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dx, double dy, double dt, double n)
	{
		psi1 = b1 * psi1 + c1 * (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx;
		psi2 = b2 * psi2 + c2 * (ex(i-1,j,k) - ex(i-1,j-1,k)) / dy;

		hz(i,j,k) -= dt / mu
				* ((ey(i,j-1,k) - ey(i-1,j-1,k)) / dx / kappa1 - (ex(i-1,j,k)
						- ex(i-1,j-1,k)) / dy / kappa2 + psi1 - psi2);
	}

protected:
	using CpmlMagnetic<T>::mu;
	using CpmlMagnetic<T>::b1;
	using CpmlMagnetic<T>::b2;
	using CpmlMagnetic<T>::c1;
	using CpmlMagnetic<T>::c2;
	using CpmlMagnetic<T>::kappa1;
	using CpmlMagnetic<T>::kappa2;
	using CpmlMagnetic<T>::psi1;
	using CpmlMagnetic<T>::psi2;
};
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_CPML_HH_*/
