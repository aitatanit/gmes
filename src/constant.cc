#include <algorithm>
#include <stdexcept>
#include "constant.hh"

using namespace std;
using namespace gmes;

const double PlusX::vector[3] = {1, 0, 0};
const double PlusY::vector[3] = {0, 1, 0};
const double PlusZ::vector[3] = {0, 0, 1};
const double MinusX::vector[3] = {-1,  0,  0};
const double MinusY::vector[3] = { 0, -1,  0};
const double MinusZ::vector[3] = { 0,  0, -1};

void PlusX::get_vector(double vector[3])
{
  copy(PlusX::vector, PlusX::vector + 3, vector);
}

void MinusX::get_vector(double vector[3])
{
  copy(MinusX::vector, MinusX::vector + 3, vector);
}

void PlusY::get_vector(double vector[3])
{
  copy(PlusY::vector, PlusY::vector + 3, vector);
}

void MinusY::get_vector(double vector[3])
{
  copy(MinusY::vector, MinusY::vector + 3, vector);
}

void PlusZ::get_vector(double vector[3])
{
  copy(PlusZ::vector, PlusZ::vector + 3, vector);
}

void MinusZ::get_vector(double vector[3])
{
  copy(MinusZ::vector, MinusZ::vector + 3, vector);
}

Component::Component()
{
  throw std::runtime_error("Component should be used as a class.");
}

Directional::Directional()
{
  throw std::runtime_error("Directional should be used as a class.");
}
