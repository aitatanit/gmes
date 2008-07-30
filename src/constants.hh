#ifndef CONSTANTS_HH_
#define CONSTANTS_HH_

#include <cmath>
#include <cfloat>


namespace gmes
{
    // physical constants

    // Vaccum permittivity
    // const double epsilon0 = 8.854187817e-12; // in SI unit (Ampere Second / Meter / Vol)
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


    class Component {};
    class Electric: public Component {};
    class Ex: public Electric {};
    class Ey: public Electric {};
    class Ez: public Electric {};
    class Magnetic: public Component {};
    class Hx: public Magnetic {};
    class Hy: public Magnetic {};
    class Hz: public Magnetic {};


    class Directional {};
    class X: public Directional {};
    class Y: public Directional {};
    class Z: public Directional {};
    class PlusX: public X {};
    class MinusX: public X {};
    class PlusY: public Y {};
    class MinusY: public Y {};
    class PlusZ: public Z {};
    class MinusZ: public Z {};
}

#endif /*CONSTANTS_HH_*/
