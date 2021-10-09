#include "Math.h"

#include <stdexcept>
#include <iostream>

namespace NumericalCalculus
{
	Point Math::Relaxation(Polynomial polynomial, double precision, Interval interval)
	{
		Polynomial firstDerivative = polynomial.takeDerivative();
		Polynomial secondDerivative = firstDerivative.takeDerivative();
		bool isTPositive = firstDerivative.isNegative(interval, precision);
		Point m = firstDerivative.findMin(interval, precision);
		Point M = firstDerivative.findMax(interval, precision);
		/*if (firstDerivative.isIncreasing(interval, precision)) {
			m = { interval.first, firstDerivative.getValue(interval.first) };
			M = { interval.second, firstDerivative.getValue(interval.second) };
		}
		else if (firstDerivative.isDecreasing(interval, precision)) {
			m = { interval.second, firstDerivative.getValue(interval.second) };
			M = { interval.first, firstDerivative.getValue(interval.first) };
		}
		else {
			m = firstDerivative.findMin(interval, precision);
			M = firstDerivative.findMax(interval, precision);
		}*/
		double t = 2 / (m.second, M.second);
		double q = (M.second - m.second) / (M.second + m.second);
		double z0 = interval.second - interval.first;
		int n = std::log(z0 / precision) / std::log(1 / q) + 1;
		double x = interval.second;
		double value = polynomial.getValue(x);
		std::cout << 0 << " " << x << " " << " " << value << std::endl;
		for (int i = 1; i <= n; i++) {
			if (isTPositive) {
				x = x + t * value;
			}
			else {
				x = x - t * value;
			}
			value = polynomial.getValue(x);
			std::cout << i << " " << x << " " << " " << value << std::endl;
		}
		return Point{ x, value };
	}
}