#include <stdexcept>

#include "pointwise_material.hh"

using namespace std;
using namespace gmes;

PointwiseMaterial::PointwiseMaterial(const int * const idx, int size)
{
    if (size != 3)
    {
    	throw domain_error("array index must be length 3.");
    }
    else
    {
	i = idx[0];
	j = idx[1];
	k = idx[2];
    }
}

