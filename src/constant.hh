#ifndef CONSTANT_HH_
#define CONSTANT_HH_

#include <cmath>

namespace gmes
{
  // physical constants

  // the ratio of the circumference of a circle to its diameter
  const double pi = 4 * atan(1);

  // the Planck constant
  const double h = 6.6260695729e-34; // in Joule/Second

  // the reduced Planck constant
  const double hbar = h / (2 * pi); // in Joule/Second

  // speed of light in vacuum
  const double c0 = 299792458; // in Meter/Second

  // vacuum permeability
  const double mu0 = 4 * pi * 1e-7; // in Henry/Meter

  // vacuum permittivity
  const double eps0 = 1 / (c0 * c0 * mu0); // in Farad/Meter

  // vacuum impedance
  const double Z0 = sqrt(mu0 / eps0); // in Ohm

  // decimal factors
  const double YOTTA = 1e24;
  const double ZETTA = 1e21;
  const double EXA   = 1e18;
  const double PETA  = 1e15;
  const double TERA  = 1e12;
  const double GIGA  = 1e9;
  const double MEGA  = 1e6;
  const double KILO  = 1e3;
  const double MILLI = 1e-3;
  const double MICRO = 1e-6;
  const double NANO  = 1e-9;
  const double PICO  = 1e-12;
  const double FEMTO = 1e-15;
  const double ATTO  = 1e-18;
  const double ZEPTO = 1e-21;
  const double YOCTO = 1e-24;

  class Component
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    
  private:
    Component();
    static const int tag = 0;
  };

  class Electric: public Component
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 1;
  };

  class Ex: public Electric
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 3;
  };

  class Ey: public Electric
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 4;
  };

  class Ez: public Electric
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 5;
  };

  class Magnetic: public Component
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 2;
  };

  class Hx: public Magnetic
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 6;
  };

  class Hy: public Magnetic
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 7;
  };

  class Hz: public Magnetic
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 8;
  };

  class ElectricCurrent: public Component
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 9;
  };

  class Jx: public ElectricCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 10;
  };

  class Jy: public ElectricCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }

  private:
    static const int tag = 11;
  };

  class Jz: public ElectricCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 12;
  };

    class MagneticCurrent: public Component
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 13;
  };

  class Mx: public MagneticCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 14;
  };

  class My: public MagneticCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 15;
  };

  class Mz: public MagneticCurrent
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 16;
  };

  class Directional
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    Directional();
    static const int tag = 17;
  };

  class X: public Directional
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 18;
  };

  class Y: public Directional
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 19;
  };

  class Z: public Directional
  {
  public:
    static int get_tag()
    {
      return tag;
    }
  private:
    static const int tag = 20;
  };

  class PlusX: public X
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 21;
  };

  class MinusX: public X
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 22;
  };

  class PlusY: public Y
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 23;
  };

  class MinusY: public Y
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 24;
  };

  class PlusZ: public Z
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 25;
  };

  class MinusZ: public Z
  {
  public:
    static int get_tag()
    {
      return tag;
    }
    static void get_vector(double vector[3]);
  private:
    static const double vector[];
    static const int tag = 26;
  };

}

#endif // CONSTANT_HH_
