#ifndef PW_DIELECTRIC_HH_
#define PW_DIELECTRIC_HH_

#include "pw_material.hh"
#include "constants.hh"

#define ex(i,j,k) ex[((this->i)*ex_y_size+(this->j))*ex_z_size+(this->k)]
#define ey(i,j,k) ey[((this->i)*ey_y_size+(this->j))*ey_z_size+(this->k)]
#define ez(i,j,k) ez[((this->i)*ez_y_size+(this->j))*ez_z_size+(this->k)]
#define hx(i,j,k) hx[((this->i)*hx_y_size+(this->j))*hx_z_size+(this->k)]
#define hy(i,j,k) hy[((this->i)*hy_y_size+(this->j))*hy_z_size+(this->k)]
#define hz(i,j,k) hz[((this->i)*hz_y_size+(this->j))*hz_z_size+(this->k)]

namespace gmes
{
template <typename T> class DielectricElectric: public MaterialElectric<T>
{
public:
	DielectricElectric(const int * const idx, int size, double epsilon = 1) :
		MaterialElectric<T>(idx, size), epsilon(epsilon)
		{
		}

	double get_epsilon()
		{
			return epsilon;
		}

	void set_epsilon(double epsilon)
		{
			epsilon = epsilon;
		}

protected:
	double epsilon;
};

template <typename T> class DielectricEx: public DielectricElectric<T>
{
public:
	DielectricEx(const int * const idx, int size, double epsilon = 1) :
			DielectricElectric<T>(idx, size, epsilon)
	{
	}

	void update(T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			double dy, double dz, double dt, double n)
	{
		ex(i,j,k) += dt / epsilon * ((hz(i+1,j+1,k) - hz(i+1,j,k)) / dy
					- (hy(i+1,j,k+1) - hy(i+1,j,k)) / dz);
	}

protected:
	using DielectricElectric<T>::epsilon;
};

template <typename T> class DielectricEy: public DielectricElectric<T>
{
public:
	DielectricEy(const int * const idx, int size, double epsilon = 1) :
			DielectricElectric<T>(idx, size, epsilon)
	{
	}

	void update(T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			double dz, double dx, double dt, double n)
	{
		ey(i,j,k) += dt / epsilon * ((hx(i,j+1,k+1) - hx(i,j+1,k)) / dz
							- (hz(i+1,j+1,k) - hz(i,j+1,k)) / dx);
	}

protected:
	using DielectricElectric<T>::epsilon;
};

template <typename T> class DielectricEz: public DielectricElectric<T>
{
public:
	DielectricEz(const int * const idx, int size, double epsilon = 1) :
			DielectricElectric<T>(idx, size, epsilon)
	{
	}

	void update(T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			double dx, double dy, double dt, double n)
	{
		ez(i,j,k) += dt / epsilon * ((hy(i+1,j,k+1) - hy(i,j,k+1)) / dx
					- (hx(i,j+1,k+1) - hx(i,j,k+1)) / dy);
	}

protected:
	using DielectricElectric<T>::epsilon;
};

template <typename T> class DielectricMagnetic: public MaterialMagnetic<T>
{
public:
	DielectricMagnetic(const int * const idx, int size, double mu = 1) :
			MaterialMagnetic<T>(idx, size), mu(mu)
		{
		}

	double get_mu()
		{
			return mu;
		}

	void set_mu(double mu)
		{
			mu = mu;
		}

protected:
	double mu;
};

template <typename T> class DielectricHx: public DielectricMagnetic<T>
{
public:
	DielectricHx(const int * const idx, int size, double mu = 1) :
			DielectricMagnetic<T>(idx, size, mu)
	{
	}

	void update(T * const hx, int hx_x_size, int hx_y_size, int hx_z_size,
			const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			double dy, double dz, double dt, double n)
	{
		hx(i,j,k) += dt / mu * ((ey(i,j-1,k) - ey(i,j-1,k-1)) / dz
				- (ez(i,j,k-1) - ez(i,j-1,k-1)) / dy);
	}

protected:
	using DielectricMagnetic<T>::mu;
};

template <typename T> class DielectricHy: public DielectricMagnetic<T>
{
public:
	DielectricHy(const int * const idx, int size, double mu = 1) :
			DielectricMagnetic<T>(idx, size, mu)
	{
	}

	void update(T * const hy, int hy_x_size, int hy_y_size, int hy_z_size,
			const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			const T * const ez, int ez_x_size, int ez_y_size, int ez_z_size,
			double dz, double dx, double dt, double n)
	{
		hy(i,j,k) += dt / mu * ((ez(i,j,k-1) - ez(i-1,j,k-1)) / dx
				- (ex(i-1,j,k) - ex(i-1,j,k-1)) / dz);
	}

protected:
	using DielectricMagnetic<T>::mu;
};

template <typename T> class DielectricHz: public DielectricMagnetic<T>
{
public:
	DielectricHz(const int * const idx, int size, double mu = 1) :
			DielectricMagnetic<T>(idx, size, mu)
	{
	}

	void update(T * const hz, int hz_x_size, int hz_y_size, int hz_z_size,
			const T * const ey, int ey_x_size, int ey_y_size, int ey_z_size,
			const T * const ex, int ex_x_size, int ex_y_size, int ex_z_size,
			double dx, double dy, double dt, double n)
	{
		hz(i,j,k) += dt / mu * ((ex(i-1,j,k) - ex(i-1,j-1,k)) / dy
				- (ey(i,j-1,k) - ey(i-1,j-1,k)) / dx);
	}

protected:
	using DielectricMagnetic<T>::mu;
};
}

#undef ex
#undef ey
#undef ez
#undef hx
#undef hy
#undef hz

#endif /*PW_DIELECTRIC_HH_*/
