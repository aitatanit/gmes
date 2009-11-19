#include <stdexcept>

#include "pw_material_real.hh"

using namespace std;
using namespace gmes;

PwMaterialReal::PwMaterialReal(const int * const idx, int size)
{
	if (size != 3)
	{
		throw domain_error("length of the array index must be 3.");
	}
	else
	{
		i = idx[0];
		j = idx[1];
		k = idx[2];
	}
}
