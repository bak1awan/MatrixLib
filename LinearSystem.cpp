#include "LinearSystem.h"
#include "Decomposition.h"

// ������� ���� ����� QR-����������
template<typename variableType>
void QRSolution(const Matrix<variableType>& A, vector<variableType>& x, const vector<variableType>& b) {
	Matrix<variableType> Q(A.rows, A.cols);
	Matrix<variableType> R(A.rows, A.cols);
	QRDecomp(A, Q, R);
	Matrix<variableType> Qt = Q.transpose();
	vector<variableType> y = Qt * b; // ������ ������ �� ������ ���� �������

	for (int i = A.cols - 1; i > -1; i--) {
		x[i] = y[i];
		for (int j = A.cols - 1; j > i; j--) {
			x[i] -= (x[j] * R[i][j]);
		}
		x[i] /= R[i][i];
	}
}

// ������� ���� ����� LU-����������
template<typename variableType>
void LUSolution(const Matrix<variableType>& A, vector<variableType>& x, const vector<variableType>& b) {
	vector<variableType> y(b.size(), 0);

	Matrix<variableType> L(A.rows, A.cols);
	Matrix<variableType> U(A.rows, A.cols);
	LUDecomp(A, L, U);

	for (int i = 0; i < A.cols; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) {
			y[i] -= (y[j] * L[i][j]);
		}
		y[i] /= L[i][i];
	}

	for (int i = A.cols - 1; i > -1; i--) {
		x[i] = y[i];
		for (int j = A.cols - 1; j > i; j--) {
			x[i] -= (x[j] * U[i][j]);
		}
		x[i] /= U[i][i];
	}
}

// ������� ���� ����� ����� ������-������
template<typename variableType>
void GaussSeidelSolution(const Matrix<variableType>& A, vector<variableType>& x, const vector<variableType>& b, variableType epsilon) {
	Matrix<variableType> B(A);
	vector<variableType> prev(x.size(), 0);

	// ��������� ������� �� ������������ ������������
	if (B.diagonal())
	{
		for (int k = 0; k < maxIteration; k++) {
			prev = x;
			for (int i = 0; i < A.rows; i++) {
				variableType var = 0;
				for (int j = 0; j < A.rows; j++)
					if (j != i) var += (B[i][j] * x[j]);

				x[i] = (b[i] - var) / B[i][i];
			}

			// �������� �� ��, �������� �� �� ������� � �������� ���������
			if (vectorLength(x - prev) < epsilon) break;

			// � ������ ������� - ������� ��������������� ���������
			if (k == (maxIteration - 1)) {
				cout << "Could not find the solution with this accuracy through " << maxIteration << " iterations.";
			}
		}

		/*
		do {
			prev = x;
			for (int i = 0; i < rows; i++) {
				double var = 0;
				for (int j = 0; j < rows; j++)
					if (j != i) var += (A[i][j] * x[j]);

				x[i] = (b[i] - var) / A[i][i];
			}
		} while (!(vectorLength(x - prev) < epsilon));
		*/

	}
	else cout << "Impossible to solve the system by this method because matrix does not have diagonal dominance.\n";
}


template void QRSolution(const Matrix<float>& A, vector<float>& x, const vector<float>& b);
template void QRSolution(const Matrix<double>& A, vector<double>& x, const vector<double>& b);
template void QRSolution(const Matrix<long double>& A, vector<long double>& x, const vector<long double>& b);


template void LUSolution(const Matrix<float>& A, vector<float>& x, const vector<float>& b);
template void LUSolution(const Matrix<double>& A, vector<double>& x, const vector<double>& b);
template void LUSolution(const Matrix<long double>& A, vector<long double>& x, const vector<long double>& b);


template void GaussSeidelSolution(const Matrix<float>& A, vector<float>& x, const vector<float>& b, float epsilon);
template void GaussSeidelSolution(const Matrix<double>& A, vector<double>& x, const vector<double>& b, double epsilon);
template void GaussSeidelSolution(const Matrix<long double>& A, vector<long double>& x, const vector<long double>& b, long double epsilon);