#pragma once

#include <map>
#include <iostream>

#include "Polynomial.h"

namespace NumericalCalculus
{
	class Math
	{
	public:
		static Point Relaxation(Polynomial polynomial, double precision, Interval interval, std::ostream& os);
	};
}
