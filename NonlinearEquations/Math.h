#pragma once

#include <map>

namespace NumericalCalculus
{
	using Polynomial = std::map<int, double>;
	class Math
	{
	public:
		static Polynomial takeDerivative(const Polynomial& polinomial);
	};
}
