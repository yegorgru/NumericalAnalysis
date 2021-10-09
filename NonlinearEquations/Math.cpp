#include "Math.h"

#include <stdexcept>

namespace NumericalCalculus
{
	double Math::Relaxation(Polynomial polynomial, double precision, Interval interval, TypeOfRoot type)
	{
		Polynomial firstDerivative = polynomial.takeDerivative();
		Polynomial secondDerivative = firstDerivative.takeDerivative();
		if (!firstDerivative.isNegative(interval, precision)) {
			throw std::runtime_error("Relaxation: first derivative is not negative");
		}
		//if(firstDerivative.isIncreasing())
	}
}