#ifndef CONSTANTS_HH_
#include <cmath>
#include <stdexcept>

namespace gmes {
// physical constants

// Vaccum permittivity
// const double epsilon0 = 8.854187817e-12; // in SI unit (Ampere Second / Meter / Vol)
#define CONSTANTS_HH_

const double epsilon0 = 1.0; // normalized

// Vacuum permeability
// const double mu0 = 1.2566370614359173e-6; // in SI unit (Second Volt / Ampere / Meter)
const double mu0 = 1.0; //  Normalized

// speed of light in vacuum
const double c0 = 1 / sqrt(epsilon0 * mu0);

// vacuum impedance
const double Z0 = sqrt(mu0 / epsilon0);

// decimal factors
const double PETA = 1e15;
const double TERA = 1e12;
const double GIGA = 1e9;
const double MEGA = 1e6;
const double KILO = 1e3;
const double MILLI = 1e-3;
const double MICRO = 1e-6;
const double NANO = 1e-9;
const double PICO = 1e-12;
const double FEMTO = 1e-15;
const double ATTO = 1e-18;

class Component {
public:
	Component() {
		throw std::runtime_error("component classes should be used as a class.");
	}
	static int get_tag() {
		return tag;
	}

protected:
	static const int tag = 0;
};

class Electric: public Component {
protected:
	static const int tag = 1;
};

class Ex: public Electric {
protected:
	static const int tag = 3;
};

class Ey: public Electric {
protected:
	static const int tag = 4;
};

class Ez: public Electric {
protected:
	static const int tag = 5;
};

class Magnetic: public Component {
protected:
	static const int tag = 2;
};

class Hx: public Magnetic {
protected:
	static const int tag = 6;
};

class Hy: public Magnetic {
protected:
	static const int tag = 7;
};

class Hz: public Magnetic {
protected:
	static const int tag = 8;
};

class Directional {
public:
	Directional() {
		throw std::runtime_error(
				"directional classes should be used as a class.");
	}
	static int get_tag() {
		return tag;
	}
	static void get_vector(double vector[3]) {
		vector[0] = vector[0];
		vector[1] = vector[1];
		vector[2] = vector[2];
	}
protected:
	static const double vector[];
	static const int tag = 10;
};

class X: public Directional {
protected:
	static const double vector[];
	static const int tag = 11;
};

class Y: public Directional {
protected:
	static const double vector[];
	static const int tag = 12;
};

class Z: public Directional {
protected:
	static const double vector[];
	static const int tag = 13;
};

class PlusX: public X {
protected:
	static const double vector[];
	static const int tag = 14;
};

class MinusX: public X {
protected:
	static const double vector[];
	static const int tag = 15;
};

class PlusY: public Y {
protected:
	static const double vector[];
	static const int tag = 16;
};

class MinusY: public Y {
protected:
	static const double vector[];
	static const int tag = 17;
};

class PlusZ: public Z {
protected:
	static const double vector[];
	static const int tag = 18;
};

class MinusZ: public Z {
protected:
	static const double vector[];
	static const int tag = 19;
};

}

#endif /*CONSTANTS_HH_*/
