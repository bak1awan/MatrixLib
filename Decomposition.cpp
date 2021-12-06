#include "Matrix.h"

// QR-разложение
template<typename variableType>
void Matrix<variableType>::QRDecomp(Matrix<variableType>& Q, Matrix<variableType>& R) {
	Matrix<variableType> AT = this->transpose();
	for (int i = 0; i < AT.cols; i++) {
		Q[i] = AT[i];
		for (int j = 0; j < i; j++)
			Q[i] = Q[i] - (scalarOp(AT[i], Q[j]) / scalarOp(Q[j], Q[j])) * Q[j];
		Q[i] = Q[i] * (1 / vectorLength(Q[i]));
	}
	R = Q * (*this);
	Q = Q.transpose();
}

// LU-разложение
template<typename variableType>
void Matrix<variableType>::LUDecomp(Matrix<variableType>& L, Matrix<variableType>& U) {
	U = *this;

	for (int k = 1; k < cols; k++)
	{
		for (int i = k - 1; i < cols; i++)
			for (int j = i; j < cols; j++)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < cols; i++)
			for (int j = k - 1; j < cols; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}
}

// Разложение Холецкого
template<typename variableType>
void Matrix<variableType>::cholesky(Matrix<variableType>& L) {
	for (int i = 0; i < rows; i++) {
		variableType res = 0;

		for (int k = 0; k < i; k++) {
			res += pow(L[i][k], 2);
		}
		L[i][i] = static_cast<variableType>(sqrt(arr[i][i] - res));

		for (int j = i + 1; j < rows; j++) {
			res = 0;

			for (int k = 0; k < i; k++) {
				res += L[i][k] * L[j][k];
			}

			L[j][i] = (arr[j][i] - res) / L[i][i];
		}
	}
}

template class Matrix<float>;
template class Matrix<double>;
template class Matrix<long double>;