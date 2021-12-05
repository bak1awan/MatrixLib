#include "Eigen.h"
#include "Decomposition.h"
#include "LinearSystem.h"

// Собственные числа через QR-разложение
template<typename variableType>
void QREigen(const Matrix<variableType>& A, vector<variableType>& x) {
	Matrix<variableType> B(A);
	Matrix<variableType> Q(A.rows, A.cols);
	Matrix<variableType> R(A.rows, A.cols);
	double error = 0;
	for (int k = 0; k < maxIteration; k++) {
		QRDecomp(B, Q, R);
		B = R * Q;
		error = 0;
		for (int i = 0; i < A.cols; i++) {
			for (int j = 0; j < A.cols; j++) {
				if (i == j) continue;
				error += pow(B[i][j], 2);
			}
		}
		if (sqrt(error) < precision) break;
	}
	for (int i = 0; i < A.cols; i++)
		x[i] = B[i][i];
}

// Собственное число через соотношение Рэлея
template<typename variableType>
variableType RayleighEigen(const Matrix<variableType>& A, variableType epsilon) {
	Matrix<variableType> B(A);
	Matrix<variableType> E(A.rows, A.cols);
	vector<variableType> x(A.rows, 1);
	vector<variableType> y(A.rows, 0);
	variableType sum{};
	variableType lambda = scalarOp(B * x, x) / scalarOp(x, x);
	variableType p{};
	for (int k = 1; k < maxIteration; k++) {
		sum = 0;
		p = lambda;
		LUSolution(B - E * lambda, y, x);
		for (int i = 0; i < A.rows; i++)
			sum += static_cast<variableType>(pow(y[i], 2));
		sum = sqrt(sum);
		for (int i = 0; i < A.rows; i++)
			x[i] = y[i] / sum;
		lambda = scalarOp(B * x, x) / scalarOp(x, x);
		if (abs(lambda - p) < epsilon) break;
	}
	return lambda;
}

// Шаблоны

template void QREigen(const Matrix<float>& A, vector<float>& x);
template void QREigen(const Matrix<double>& A, vector<double>& x);
template void QREigen(const Matrix<long double>& A, vector<long double>& x);

template float RayleighEigen(const Matrix<float>& A, float epsilon);
template double RayleighEigen(const Matrix<double>& A, double epsilon);
template long double RayleighEigen(const Matrix<long double>& A, long double epsilon);
