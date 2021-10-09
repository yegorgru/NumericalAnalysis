#include "Math.h"

namespace NumericalCalculus
{
	Polynomial Math::takeDerivative(const Polynomial& polinomial)
	{
		Polynomial derivative;
		for (const auto& [degree, coef] : polinomial) {
			if (degree != 0) {
				derivative[degree - 1] = degree * coef;
			}
		}
		return derivative;
	}
}