#include "LinearSystemSolver.h"

#include <stdexcept>
#include <cmath>

namespace NumericalAnalysis {

	std::vector<double> LinearSystemSolver::GaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b)
	{
		if (A.size() == 0 || A.at(0).size() == 0) {
			return {};
		}
		CheckRange(A);
		CheckIndependence(A);
		for (size_t cur = 0; cur < A.size(); ++cur) {
			for (size_t i = cur; i < A.size(); ++i) {
				if (A.at(i).at(cur) != 0.0) {
					A.at(cur).swap(A.at(i));
					std::swap(b.at(cur), b.at(i));
					break;
				}
			}
			double coef = A.at(cur).at(cur);
			for (size_t i = 0; i < A.at(cur).size(); ++i) {
				A.at(cur).at(i) /= coef;
			}
			b.at(cur) /= coef;
			for (size_t i = 0; i < A.at(cur).size(); ++i) {
				if (i == cur || A.at(i).at(cur) == 0) {
					continue;
				}
				double coef = A.at(i).at(cur);
				for (size_t k = 0; k < A.at(i).size(); ++k) {
					A.at(i).at(k) -= coef * A.at(cur).at(k);
				}
				b.at(i) -= coef * b.at(cur);
			}
		}
		return b;
	}

	void LinearSystemSolver::CheckRange(const std::vector<std::vector<double>>& A)
	{
		if (A.at(0).size() != A.size()) {
			throw std::runtime_error("Incorrect size of A");
		}
		if (A.size() == 0) {
			return;
		}
		size_t size = A.at(0).size();
		for (size_t i = 1; i < A.size(); ++i) {
			if (A.at(i).size() != size) {
				throw std::runtime_error("Different length of vectors");
			}

		}
	}

	void LinearSystemSolver::CheckIndependence(const std::vector<std::vector<double>>& A)
	{
		if (A.size() == 0 || A.at(0).size() == 0) {
			return;
		}
		for (size_t i = 0; i < A.size(); ++i) {
			for (size_t j = i + 1; j < A.size(); ++j) {
				double coef = INT_MAX;
				for (size_t k = 0; k < A.at(i).size(); ++k) {
					if (A.at(j).at(k) == 0) {
						if (A.at(i).at(k) == 0) {
							continue;
						}
						else {
							break;
						}
					}
					if (coef == INT_MAX) {
						coef = A.at(i).at(k) / A.at(j).at(k);
						continue;
					}
					if (A.at(i).at(k) / A.at(j).at(k) != coef) {
						break;
					}
					if (k == A.at(i).size() - 1) {
						throw std::runtime_error("Dependent vectors");
					}
				}
			}
		}
	}
}