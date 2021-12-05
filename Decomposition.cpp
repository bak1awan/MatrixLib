#include "Decomposition.h"

// QR-разложение
template<typename variableType>
void QRDecomp(const Matrix<variableType>& A, Matrix<variableType>& Q, Matrix<variableType>& R) {
	Matrix<variableType> AT = A.transpose();
	for (int i = 0; i < AT.cols; i++) {
		Q[i] = AT[i];
		for (int j = 0; j < i; j++)
			Q[i] = Q[i] - (scalarOp(AT[i], Q[j]) / scalarOp(Q[j], Q[j])) * Q[j];
		Q[i] = Q[i] * (1 / vectorLength(Q[i]));
	}
	R = Q * A;
	Q = Q.transpose();
}

// LU-разложение
template<typename variableType>
void LUDecomp(const Matrix<variableType>& A, Matrix<variableType>& L, Matrix<variableType>& U) {
	U = A;

	for (int k = 1; k < A.cols; k++)
	{
		for (int i = k - 1; i < A.cols; i++)
			for (int j = i; j < A.cols; j++)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < A.cols; i++)
			for (int j = k - 1; j < A.cols; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}
}

// Разложение Холецкого
template<typename variableType>
void cholesky(const Matrix<variableType>& A, Matrix<variableType>& L) {
	for (int i = 0; i < A.rows; i++) {
		variableType res = 0;

		for (int k = 0; k < i; k++) {
			res += pow(L[i][k], 2);
		}
		L[i][i] = static_cast<variableType>(sqrt(A[i][i] - res));

		for (int j = i + 1; j < A.rows; j++) {
			res = 0;

			for (int k = 0; k < i; k++) {
				res += L[i][k] * L[j][k];
			}

			L[j][i] = (A[j][i] - res) / L[i][i];
		}
	}
}


// Шаблоны

template void QRDecomp(const Matrix<float>& A, Matrix<float>& Q, Matrix<float>& R);
template void QRDecomp(const Matrix<double>& A, Matrix<double>& Q, Matrix<double>& R);
template void QRDecomp(const Matrix<long double>& A, Matrix<long double>& Q, Matrix<long double>& R);

template void LUDecomp(const Matrix<float>& A, Matrix<float>& L, Matrix<float>& U);
template void LUDecomp(const Matrix<double>& A, Matrix<double>& L, Matrix<double>& U);
template void LUDecomp(const Matrix<long double>& A, Matrix<long double>& L, Matrix<long double>& U);

template void cholesky(const Matrix<float>& A, Matrix<float>& L);
template void cholesky(const Matrix<double>& A, Matrix<double>& L);
template void cholesky(const Matrix<long double>& A, Matrix<long double>& L);
