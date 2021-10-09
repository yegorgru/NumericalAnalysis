#include "Math.h"

#include <stdexcept>

namespace NumericalCalculus
{
	Iterations Math::Relaxation(Polynomial polynomial, double precision, Interval interval)
	{
		Polynomial firstDerivative = polynomial.takeDerivative();
		bool isTPositive = firstDerivative.isNegative(interval, precision);
		Point m = firstDerivative.findAbsMin(interval, precision);
		Point M = firstDerivative.findAbsMax(interval, precision);
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
		double t = 2 / (m.second + M.second);
		double q = (M.second - m.second) / (M.second + m.second);
		double z0 = interval.second - interval.first;
		int n = std::log(z0 / precision) / std::log(1 / q) + 1;
		double x = interval.second;
		double value = polynomial.getValue(x);
		Iterations result{ {x, value} };
		for (int i = 0; i < n; i++) {
			if (isTPositive) {
				x = x + t * value;
			}
			else {
				x = x - t * value;
			}
			value = polynomial.getValue(x);
			result.emplace_back(x, value);
		}
		return result;
	}

	Iterations Math::Newton(Polynomial polynomial, double precision, Interval interval)
	{
		if (polynomial.getValue(interval.first) * polynomial.getValue(interval.second) >= 0) {
			throw std::runtime_error("values on sides of interval are not with different signs");
		}
		Polynomial firstDerivative = polynomial.takeDerivative();
		Polynomial secondDerivative = firstDerivative.takeDerivative();
		if (secondDerivative.changeSign(interval, precision)) {
			throw std::runtime_error("second derivative change sign");
		}
		Point m = firstDerivative.findAbsMin(interval, precision);
		Point M = secondDerivative.findAbsMax(interval, precision);
		double x = interval.second;
		double maxMistake = std::max(std::abs(x - interval.second), std::abs(x - interval.first));
		double q = M.second * maxMistake / 2 / m.second;
		while (polynomial.getValue(x) * secondDerivative.getValue(x) <= 0 || q >= 1) {
			x -= precision;
			maxMistake = std::max(std::abs(x - interval.second), std::abs(x - interval.first));
			q = M.second * maxMistake / 2 / m.second;
			if (x < interval.first) {
				throw std::runtime_error("the conditions of the theorem cannot be satisfied");
			}
		}
		int n = std::log2(
							(std::log((interval.second - interval.first) / precision) 
								/ 
							std::log(1 / q) + 1)
						) + 1;
		double value = polynomial.getValue(x);
		Iterations result{ {x, value} };
		for (int i = 0; i < n; i++) {
			x = x - value / firstDerivative.getValue(x);
			value = polynomial.getValue(x);
			result.emplace_back(x, value);
		}
		return result;
	}

	Iterations Math::Secant(Polynomial polynomial, double precision, Interval interval)
	{
		if (polynomial.getValue(interval.first) * polynomial.getValue(interval.second) >= 0) {
			throw std::runtime_error("values on sides of interval are not with different signs");
		}
		Polynomial firstDerivative = polynomial.takeDerivative();
		Polynomial secondDerivative = firstDerivative.takeDerivative();
		if (secondDerivative.changeSign(interval, precision)) {
			throw std::runtime_error("second derivative change sign");
		}
		Point m = firstDerivative.findAbsMin(interval, precision);
		Point M = secondDerivative.findAbsMax(interval, precision);
		double x = interval.second;
		double maxMistake = std::max(std::abs(x - interval.second), std::abs(x - interval.first));
		double q = M.second * maxMistake / 2 / m.second;
		while (polynomial.getValue(x) * secondDerivative.getValue(x) <= 0 || q >= 1) {
			x -= precision;
			maxMistake = std::max(std::abs(x - interval.second), std::abs(x - interval.first));
			q = M.second * maxMistake / 2 / m.second;
			if (x < interval.first) {
				throw std::runtime_error("the conditions of the theorem cannot be satisfied");
			}
		}
		double x1 = x - precision;
		maxMistake = std::max(std::abs(x1 - interval.second), std::abs(x1 - interval.first));
		q = M.second * maxMistake / 2 / m.second;
		while (polynomial.getValue(x1) * secondDerivative.getValue(x1) <= 0 || q >= 1) {
			x1 -= precision;
			maxMistake = std::max(std::abs(x1 - interval.second), std::abs(x1 - interval.first));
			q = M.second * maxMistake / 2 / m.second;
			if (x1 < interval.first) {
				throw std::runtime_error("the conditions of the theorem cannot be satisfied");
			}
		}
		int n = std::log2(
			(std::log((interval.second - interval.first) / precision)
				/
				std::log(1 / q) + 1)
		) + 1;
		double value = polynomial.getValue(x);
		double value1 = polynomial.getValue(x1);
		Iterations result{ {x, value}, {x1, value1} };
		for (int i = 0; i < n; i++) {
			double nextX = x1 - (x1 - x) * value1
								/
								(value1 - value);
			double nextValue = polynomial.getValue(nextX);
			x = x1;
			value = value1;
			x1 = nextX;
			value1 = nextValue;
			result.emplace_back(x1, value1);
			if (std::abs(value1) < precision) {
				result.resize(n, { x, value1 });
				break;
			}
		}
		return result;
	}
}