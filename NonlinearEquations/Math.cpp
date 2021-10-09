#include "Math.h"

namespace NumericalCalculus
{
	double Math::Relaxation(Polynomial polynomial, double precision, Interval interval, TypeOfRoot type)
	{
		Polynomial firstDerivative = polynomial.takeDerivative();
		Polynomial secondDerivative = firstDerivative.takeDerivative();
		return 0.0;
	}
}