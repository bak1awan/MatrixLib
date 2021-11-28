#include "Matrix.h"
#include <math.h>
#include <iomanip>

const double precision = 1e-5;
const int maxIteration = 40;

Matrix::Matrix(int nRow, int nCol) : rows(nRow), cols(nCol) {
	arr = vector<vector<double>>(nRow, vector <double >(nCol, 0));
}

// ������� �����
Matrix Matrix::operator- () 
{
	Matrix result(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result[i][j] = -result[i][j];
		}
	}
	return result;
}

// ���������������� �������
Matrix Matrix::transpose() {
	Matrix result(cols, rows);
	for (int i = 0; i < result.cols; i++) {
		for (int j = 0; j < result.rows; j++)
			result[i][j] = arr[j][i];
	}
	return result;
}

// �������� ������
Matrix Matrix::operator+ (const Matrix& rhs) {
	Matrix result(*this);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] += rhs[i][j];
		}
	}
	return result;
}

// ��������� ������
Matrix Matrix::operator- (const Matrix& rhs) {
	Matrix result(*this);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] -= rhs[i][j];
		}
	}
	return result;
}

// ��������� ������
Matrix Matrix::operator* (const Matrix& rhs) {
	Matrix result(rows, rhs.cols);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			for (int r = 0; r < cols; r++)
				result[i][j] += arr[i][r] * rhs[r][j];
			if (abs(result[i][j]) < precision)
				result[i][j] = 0;
		}
	}
	return result;
}

// ��������� ������� �� ������ ������
vector<double> operator* (Matrix& A, const vector<double>& x) {
	vector<double> result(x.size(), 0);
	for (int i = 0; i < A.cols; i++) {
		for (int j = 0; j < A.cols; j++)
			result[i] += A[i][j] * x[j];
	}
	return result;
}

// ��������� ������� �� ������ �����
vector<double> operator* (const vector<double>& x, Matrix& A) {
	return A * x;
}

// ��������� ������� �� ����� ������
Matrix operator* (Matrix &A, double n) {
	Matrix result(A);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] *= n;
		}
	}
	return result;
}

// ��������� ������� �� ����� �����
Matrix operator* (double n, Matrix& A) {
	Matrix result(A);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] *= n;
		}
	}
	return result;
}

// ���������� ��������� ������
std::ostream& operator << (std::ostream& out, const Matrix& A) {
	for (int i = 0; i < A.arr.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			out << A[i][j] << "\t";
		}
		out << endl;
	}
	return out;
}

// ���������� ��������� �����
std::istream& operator >> (std::istream& in, Matrix& A) {
	for (int i = 0; i < A.arr.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			cout << " Input [" << i << "]" << "["
				<< j << "]" << " element : ";
			in >> A[i][j];
		}
	}
	return in;
}

// ��������� ������������ ��������
double scalarOp(const vector<double>& v1, const vector<double>& v2) {
	double result(0);
	for (int i = 0; i < v1.size(); i++)
		result += v1[i] * v2[i];
	return result;
}

// QR-����������
void Matrix::QRDecomp(Matrix& Q, Matrix& R) {
	Matrix AT = this->transpose();
	for (int i = 0; i < AT.cols; i++) {
		Q[i] = AT[i];
		for (int j = 0; j < i; j++)
			Q[i] = Q[i] - (scalarOp(AT[i], Q[j]) / scalarOp(Q[j], Q[j])) * Q[j];
		Q[i] = Q[i] * (1 / vectorLength(Q[i]));
	}
	R = Q * (*this);
	Q = Q.transpose();
}

// LU-����������
void Matrix::LUDecomp(Matrix& L, Matrix& U) {
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

// ������������ ����� LU-����������
double Matrix::LUDeterminant() {
	double result = 1;

	Matrix L(rows, cols);
	Matrix U(rows, cols);

	this->LUDecomp(L, U);

	for (int i = 0; i < cols; i++)
		result *= U[i][i];
	return result;
}

// �������� ���������� ��� ������
vector<double>& Matrix::operator[] (const int i) {
	return arr[i];
}

const vector<double>& Matrix::operator[] (const int i) const {
	return arr[i];
}

void Matrix::cholesky(Matrix& L) {
	for (int i = 0; i < rows; i++) {
		double res = 0;

		for (int k = 0; k < i; k++) {
			res += pow(L[i][k], 2);
		}

		L[i][i] = sqrt(arr[i][i] - res);

		for (int j = i + 1; j < rows; j++) {
			res = 0;

			for (int k = 0; k < i; k++) {
				res += L[i][k] * L[j][k];
			}

			L[j][i] = (arr[j][i] - res) / L[i][i];
		}
	}
}

// ��������� ������� �� ����� �����
vector<double> operator* (double n, const vector<double>& v) {
	vector<double> result(v.size(), 0);
	for (int i = 0; i < v.size(); i++) {
		if (abs(n * v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n * v[i];
	}
	return result;
}

// ��������� ������� �� ����� ������
vector<double> operator* (const vector<double>& v, double n) {
	return n * v;
}

// ��������� ����� �� �������
vector<double> operator- (const vector<double>& v, double n) {
	vector<double> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] - n;
	}
	return result;
}


// �������� ������� � ������ �����
vector<double> operator+ (double n, const vector<double>& v) {
	vector<double> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] + n;
	}
	return result;
}

// �������� ������� � ������ ������
vector<double> operator+ (const vector<double>& v, double n) {
	return n + v;
}

// �������� ��������
vector<double> operator+ (const vector<double>& v1, const vector<double>& v2) {
	vector<double> result(v1.size(), 0);
	for (int i = 0; i < v1.size(); i++)
		result[i] = v1[i] + v2[i];
	return result;
}

// ��������� ��������
vector<double> operator- (const vector<double>& v1, const vector<double>& v2) {
	vector<double> result(v1.size(), 0);
	for (int i = 0; i < v1.size(); i++)
		result[i] = v1[i] - v2[i];
	return result;
}

// ������� ���� ����� LU-����������
void Matrix::LUSolution(vector<double>& x, const vector<double>& b) {
	vector<double> y(b.size(), 0);

	Matrix L(rows, cols);
	Matrix U(rows, cols);
	this->LUDecomp(L, U);

	for (int i = 0; i < cols; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) {
			y[i] -= (y[j] * L[i][j]);
		}
		y[i] /= L[i][i];
	}

	for (int i = cols - 1; i > -1; i--) {
		x[i] = y[i];
		for (int j = cols - 1; j > i; j--) {
			x[i] -= (x[j] * U[i][j]);
		}
		x[i] /= U[i][i];
	}
}

// ������� ���� ����� QR-����������
void Matrix::QRSolution(vector<double>& x, const vector<double>& b) {
	Matrix Q(rows, cols);
	Matrix R(rows, cols);
	this->QRDecomp(Q, R);
	Matrix Qt = Q.transpose();
	vector<double> y = Qt * b; // ������ ������ �� ������ ���� �������


	for (int i = cols - 1; i > -1; i--) {
		x[i] = y[i];
		for (int j = cols - 1; j > i; j--) {
			x[i] -= (x[j] * R[i][j]);
		}
		x[i] /= R[i][i];
	}
}

// ����������� ����� ����� QR-����������
void Matrix::QREigen(vector<double>& x) {
	Matrix B(rows, cols);
	B = *this;
	Matrix Q(rows, cols);
	Matrix R(rows, cols);
	double error = 0;
	for (int k = 0; k < maxIteration; k++) {
		B.QRDecomp(Q, R);
		B = R * Q;
		error = 0;
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < cols; j++) {
				if (i == j) continue;
				error += pow(B[i][j], 2);
			}
		}
		if (sqrt(error) < precision) break;
	}
	for (int i = 0; i < cols; i++)
		x[i] = B[i][i];
}

// ���������� ����� ������� ����� ��� ����������
double vectorLength(const vector<double>& v) {
	double length = 0;
	for (int i = 0; i < v.size(); i++)
		length += std::pow(v[i], 2);
	return std::sqrt(length);
}

// ��������� ��������� �������
Matrix generateE(int size) {
	Matrix E(size, size);
	for (int i = 0; i < size; i++)
		E[i][i] = 1;
	return E;
}

// �������� ������� ����� LU-����������
Matrix Matrix::LUInverse() {
	// ���������� ��������� ������� 
	Matrix E = generateE(rows);

	// ������� ��� �����������
	Matrix X(rows, cols);

	// � ����� �������� ������� ���� ����� LU-���������� �� �������� ������ X � E
	for (int i = 0; i < cols; i++)
		this->LUSolution(X[i], E[i]);

	// ������������� ������� X, ����� �������� ���������� ��� �������
	return X.transpose();
}

// �������� ������� ����� QR-����������
Matrix Matrix::QRInverse() {
	// ���������� ������� ��� ����������
	Matrix Q(rows, cols);
	Matrix R(rows, cols);

	// �������� ������� Q � R
	this->QRDecomp(Q, R);

	// ������� �����������
	Matrix X(rows, cols);

	// ������� ��� �������
	vector<double> y(rows, 0);
	vector<double> x(rows, 0);

	// �������� ������ ��� ������� X, ��� ��� ��������� ���� ������ �� �������
	for (int k = 0; k < cols; k++) {
		x = X[k];
		y = Q[k];
		for (int i = cols - 1; i > -1; i--) {
			x[i] = y[i];
			for (int j = cols - 1; j > i; j--) {
				x[i] -= (x[j] * R[i][j]);
			}
			x[i] /= R[i][i];
		}
		X[k] = x;
	}

	// ������������� ���������� �������, ����� �������� ���������� ���
	return X.transpose();
}

// ����� �������
double Matrix::norm() {
	double result = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++)
			result += pow((*this)[i][j], 2);
	}
	return sqrt(result);
}

// �������� ������� ������� ������ (Newton-Schulz-Hotelling)
Matrix Matrix::ShultzInverse(double epsilon) {

	Matrix A = *this;
	Matrix E = generateE(rows);
	Matrix phi(rows, cols);
	Matrix At = A.transpose();
	// �������� ��������� ������� X
	Matrix U = At * (1.0 / (A * At).norm()); // ����� ���� ������� ��������

	for (int i = 0;; i++) {
		phi = E - A * U;
		if (phi.norm() < epsilon) break;
		U = U * (E + phi);
	}

	return U;
}

// ������� ���� ����� ����� ������-������
void Matrix::GaussSeidelSolution(vector<double>& x, const vector<double>& b, double epsilon) {
	Matrix A = *this;
	vector<double> prev(x.size(), 0);

	// ��������� ������� �� ������������ ������������
	if (A.diagonal())
	{
		for (int k = 0; k < maxIteration; k++) {
			prev = x;
			for (int i = 0; i < rows; i++) {
				double var = 0;
				for (int j = 0; j < rows; j++)
					if (j != i) var += (A[i][j] * x[j]);

				x[i] = (b[i] - var) / A[i][i];
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

// �������� ������� �� ������������ ������������
bool Matrix::diagonal() {
	Matrix A = *this;
	double sum;
	for (int i = 0; i < rows; i++) {
		sum = 0;
		for (int j = 0; j < rows; j++) 
			sum += abs(A[i][j]);

		sum -= abs(A[i][i]);

		if (sum > A[i][i]) return false;
	}
	return true;
}

// ����������� ����� ����� ����������� �����
double Matrix::RayleighEigen(double epsilon) {
	Matrix A = *this;
	Matrix E = generateE(rows);
	vector<double> x(rows, 1);
	vector<double> y(rows, 0);
	double sum;
	double lambda = scalarOp(A * x, x) / scalarOp(x, x);
	double p;
	for (int k = 1; k < maxIteration; k++) {
		sum = 0;
		p = lambda;
		(A - E * lambda).LUSolution(y, x);
		for (int i = 0; i < rows; i++)
			sum += pow(y[i], 2);
		sum = sqrt(sum);
		for (int i = 0; i < rows; i++)
			x[i] = y[i] / sum;
		lambda = scalarOp(A * x, x) / scalarOp(x, x);
		if (abs(lambda - p) < epsilon) break;
	}
	return lambda;
}

// �������� ������� � ������ ������
Matrix operator+ (Matrix& A, double x) {
	Matrix result(A.rows, A.cols);
	for (int i = 0; i < A.rows; i++)
		for (int j = 0; j < A.cols; j++)
			result[i][j] = A[i][j] + x;
	return result;
}

// �������� ������� � ������ �����
Matrix operator+ (double x, Matrix& A) {
	return A + x;
}

Matrix operator* (const vector<double>& v1, const vector<double>& v2) {
	Matrix res(v1.size(), v1.size());
	for (int i = 0; i < v1.size(); i++) {
		for (int j = 0; j < v1.size(); j++)
			res[i][j] = v1[i] * v2[j];
	}
	return res;
}

double& VectorT::operator[](int i) {
	return v[i];
}

const double& VectorT::operator[](int i) const {
	return v[i];
}

int VectorT::size() {
	return this->v.size();
}

// ��������� ������������������ ������� �� ������� - �����
double VectorT::operator* (vector<double>& v1) {
	double result(0);
	for (int i = 0; i < v.size(); i++)
		result += v[i] * v1[i];
	return result;
}

// ����� ����������������� ��������
VectorT VectorT::operator+ (VectorT& v1) {
	VectorT result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] + v1[i];
	return result;
}

// �������� ����������������� ��������
VectorT VectorT::operator- (VectorT& v1) {
	VectorT result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] - v1[i];
	return result;
}

// �������� � ������
VectorT  VectorT::operator+ (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] + n;
	}
	return result;
}

VectorT operator+ (double n, VectorT& v) {
	return v + n;
}

// ��������� �����
VectorT VectorT::operator- (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] - n;
	}
	return result;
}

// ��������� �� �����
VectorT VectorT::operator* (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n * v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n * v[i];
	}
	return result;
}

VectorT operator* (double n, VectorT& v) {
	return v * n;
}

// ������� �� �����
VectorT VectorT::operator/ (double n) {
	VectorT result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n / v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n / v[i];
	}
	return result;
}

// ����� ������������������ �������
ostream& operator<< (ostream& out , VectorT& vT) {
	for (double el : vT.v) {
		out << el << ' ';
	}
	return out;
}

// ���������� ������������������ �������
istream& operator>> (istream& in, VectorT& vT) {
	for (int i = 0; i < vT.size(); i++) {
		cout << " Input [" << i << "]" << " element : ";
		in >> vT[i];
	}
	return in;
}

Matrix operator* (const vector<double>& v1, const VectorT& v2) {
	Matrix res(v1.size(), v1.size());
	for (int i = 0; i < v1.size(); i++) {
		for (int j = 0; j < v1.size(); j++)
			res[i][j] = v1[i] * v2[j];
	}
	return res;
}

VectorT transpose(vector<double>& v) {
	return VectorT(v);
}