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

	std::vector<double> LinearSystemSolver::TridiagonalMatrix(std::vector<std::vector<double>>& A, std::vector<double>& b)
	{
		if (A.size() == 0 || A.at(0).size() == 0) {
			return {};
		}
		CheckRange(A);
		CheckIndependence(A);
		CheckTridiagonal(A);

		for (size_t i = 1; i < A.size(); ++i) {
			double w = A.at(i).at(i - 1) / A.at(i - 1).at(i - 1);
			A.at(i).at(i) -= w * A.at(i - 1).at(i);
			b.at(i) -= w * b.at(i - 1);
		}
		std::vector<double> answer(A.size(), 0);
		answer.back() = b.back() / A.back().back();
		for (int i = answer.size() - 2; i >= 0; --i) {
			answer.at(i) = (b.at(i) - A.at(i).at(i + 1) * answer.at(i + 1)) / A.at(i).at(i);
		}
		return answer;
	}

	std::vector<std::vector<double>> LinearSystemSolver::Jacobi(std::vector<std::vector<double>>& A, std::vector<double>& b, double precision)
	{
		if (A.size() == 0 || A.at(0).size() == 0) {
			return {};
		}
		CheckRange(A);
		CheckIndependence(A);
		CheckZeroOnMainDiagonal(A);

		double q = -1.0;
		for (size_t i = 0; i < A.size(); ++i) {
			double sum = 0;
			for (size_t j = 0; j < A.at(i).size(); ++j) {
				if (i != j) {
					sum += std::abs(A.at(i).at(j));
				}
			}
			q = std::max(q, sum / std::abs(A.at(i).at(i)));
		}
		/*if (q >= 1.0) {
			throw std::runtime_error("q >= 1, not convergent");
		}*/
		int n = q >= 1.0 ? 30 : std::log(precision * (1 - q)) / std::log(q) + 1;
		
		std::vector<double> answer(b.size(), 0);
		std::vector<std::vector<double>> result{answer};
		std::vector<double> nextAnswer(b.size(), 0);
		for (int counter = 0; counter < n; ++counter) {
			for (size_t i = 0; i < nextAnswer.size(); ++i) {
				double sum = 0;
				for (size_t j = 0; j < nextAnswer.size(); ++j) {
					if (i != j) {
						sum += A.at(i).at(j) * answer.at(j);
					}
				}
				nextAnswer.at(i) = 1 / A.at(i).at(i) * (b.at(i) - sum);
			}
			answer = nextAnswer;
			result.push_back(nextAnswer);
		}
		return result;
	}

	std::vector<std::vector<double>> LinearSystemSolver::UpperRelaxation(std::vector<std::vector<double>>& A, std::vector<double>& b, double precision)
	{
		if (A.size() == 0 || A.at(0).size() == 0) {
			return {};
		}
		CheckRange(A);
		CheckIndependence(A);
		CheckZeroOnMainDiagonal(A);

		double q = -1.0;
		for (size_t i = 0; i < A.size(); ++i) {
			double sum = 0;
			for (size_t j = 0; j < A.at(i).size(); ++j) {
				if (i != j) {
					sum += std::abs(A.at(i).at(j));
				}
			}
			q = std::max(q, sum / std::abs(A.at(i).at(i)));
		}
		/*if (q >= 1.0) {
			throw std::runtime_error("q >= 1, not convergent");
		}*/
		int n = 20;
		std::vector<double> answer(b.size(), 0);
		std::vector<std::vector<double>> result{ answer };
		std::vector<double> nextAnswer(b.size(), 0);
		double w = 1.5;
		for (int counter = 0; counter < n; ++counter) {
			for (size_t i = 0; i < nextAnswer.size(); ++i) {
				double sum = 0;
				for (size_t j = 0; j < i; ++j) {
					sum += A.at(i).at(j) * nextAnswer.at(j);
				}
				for (size_t j = i + 1; j < nextAnswer.size(); ++j) {
					sum += A.at(i).at(j) * answer.at(j);
				}
				nextAnswer.at(i) = (1 - w) * answer.at(i) + 1 / A.at(i).at(i) * w * (b.at(i) - sum);
			}
			answer = nextAnswer;
			result.push_back(nextAnswer);
		}
		return result;
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

	void LinearSystemSolver::CheckTridiagonal(const std::vector<std::vector<double>>& A)
	{
		if (A.size() < 3) {
			return;
		}
		for (size_t i = 0; i < A.size(); ++i) {
			for (size_t j = 0; j < A.at(i).size(); ++j) {
				if (std::max(i, j) - std::min(i, j) > 1 && A.at(i).at(j) != 0.0) {
					throw std::runtime_error("A is not tridiagonal");
				}
			}
		}
	}

	void LinearSystemSolver::CheckZeroOnMainDiagonal(const std::vector<std::vector<double>>& A) {
		for (size_t i = 0; i < A.size(); ++i) {
			if (A.at(i).at(i) == 0) {
				throw std::runtime_error("main diagonal has 0");
			}
		}
	}

}

