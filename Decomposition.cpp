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

// LDLT-разложение
template<typename variableType>
void Matrix<variableType>::LDLDecomp(Matrix<variableType>& L, Matrix<variableType>& D) {
	// Копируем матрицу во внутреннюю переменную
	Matrix<variableType> A(*this);

	// Алгоритм для получения элементов матриц L и D
	for (int k = 0; k < rows - 1; k++) {
		for (int i = k + 1; i < rows; i++) {
			A[i][k] /= A[k][k];
			for (int j = k + 1; j < rows; j++) {
				A[i][j] -= A[i][k] * A[k][j];
			}
		}
	}

	// Записываем по главной диагонали матрицы D элементы с главной диагонали полученной матрицы A
	// а матрице L по главной диагонали выставляем единицы
	for (int i = 0; i < rows; i++) {
		D[i][i] = A[i][i];
		L[i][i] = 1;
	}

	// Заполняем элементы ниже главной диагонали матрицы L элементами ниже главной диагонали из матрицы A
	for (int i = rows / 2; i < rows; i++) {
		for (int j = 0; j < i; j++) {
			L[i][j] = A[i][j];
		}
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