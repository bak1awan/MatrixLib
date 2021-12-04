#include "Matrix.h"
#include <math.h>
#include <iomanip>

const double precision = 1e-5;
const int maxIteration = 40;

// ������� �����
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator- ()
{
	Matrix<variableType> result(rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result[i][j] = -result[i][j];
		}
	}
	return result;
}

// ���������������� �������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::transpose() {
	Matrix<variableType> result(cols, rows);
	for (int i = 0; i < result.cols; i++) {
		for (int j = 0; j < result.rows; j++)
			result[i][j] = arr[j][i];
	}
	return result;
}

// �������� ������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator+ (const Matrix<variableType>& rhs) {
	Matrix<variableType> result(*this);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] += rhs[i][j];
		}
	}
	return result;
}

// ��������� ������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator- (const Matrix<variableType>& rhs) {
	Matrix<variableType> result(*this);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] -= rhs[i][j];
		}
	}
	return result;
}

// ��������� ������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator* (const Matrix<variableType>& rhs) {
	Matrix<variableType> result(rows, cols);
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
template<typename variableType>
vector<variableType> Matrix<variableType>::operator* (const vector<variableType>& x) {
	vector<variableType> result(x.size(), 0);
	for (int i = 0; i < cols; i++) {
		for (int j = 0; j < cols; j++)
			result[i] += arr[i][j] * x[j];
	}
	return result;
}

// ��������� ������� �� ������ �����
template<typename variableType>
vector<variableType> operator* (const vector<variableType>& x, Matrix<variableType>& A) {
	return A * x;
}

// ��������� ������� �� ����� ������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator* (variableType n) {
	Matrix<variableType> result(*this);
	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			result[i][j] *= n;
		}
	}
	return result;
}

// ��������� ������� �� ����� �����
template<typename variableType>
Matrix<variableType> operator* (variableType n, Matrix<variableType>& A) {
	return A * n;
}

// �������� ������� � ������ ������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::operator+ (variableType x) {
	Matrix<variableType> result(rows, cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result[i][j] = arr[i][j] + x;
	return result;
}

// �������� ������� � ������ �����
template<typename variableType>
Matrix<variableType> operator+ (variableType x, Matrix<variableType>& A) {
	return A + x;
}

// ���������� ��������� ������ ��� ������
template<typename variableType>
std::ostream& operator << (std::ostream& out, const Matrix<variableType>& A) {
	for (int i = 0; i < A.arr.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			out << A[i][j] << "\t";
		}
		out << endl;
	}
	return out;
}

// ���������� ��������� ����� ��� ������
template<typename variableType>
std::istream& operator >> (std::istream& in, Matrix<variableType>& A) {
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
template<typename variableType>
variableType scalarOp(const vector<variableType>& v1, const vector<variableType>& v2) {
	variableType result(0);
	for (int i = 0; i < v1.size(); i++)
		result += v1[i] * v2[i];
	return result;
}

// QR-����������
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

// LU-����������
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

// ������������ ����� LU-����������
template<typename variableType>
variableType Matrix<variableType>::LUDeterminant() {
	variableType result = 1;

	Matrix<variableType> L(rows, cols);
	Matrix<variableType> U(rows, cols);

	this->LUDecomp(L, U);

	for (int i = 0; i < cols; i++)
		result *= U[i][i];
	return result;
}

// �������� ���������� ��� ������
template<typename variableType>
vector<variableType>& Matrix<variableType>::operator[] (const int i) {
	return arr[i];
}

// �������� ���������� ��� ����������� ������
template<typename variableType>
const vector<variableType>& Matrix<variableType>::operator[] (const int i) const {
	return arr[i];
}

// ���������� ���������
template<typename variableType>
void Matrix<variableType>::cholesky(Matrix<variableType>& L) {
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

// ������� ���� ����� LU-����������
template<typename variableType>
void Matrix<variableType>::LUSolution(vector<variableType>& x, const vector<variableType>& b) {
	vector<variableType> y(b.size(), 0);

	Matrix<variableType> L(rows, cols);
	Matrix<variableType> U(rows, cols);
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
template<typename variableType>
void Matrix<variableType>::QRSolution(vector<variableType>& x, const vector<variableType>& b) {
	Matrix<variableType> Q(rows, cols);
	Matrix<variableType> R(rows, cols);
	this->QRDecomp(Q, R);
	Matrix<variableType> Qt = Q.transpose();
	vector<variableType> y = Qt * b; // ������ ������ �� ������ ���� �������

	for (int i = cols - 1; i > -1; i--) {
		x[i] = y[i];
		for (int j = cols - 1; j > i; j--) {
			x[i] -= (x[j] * R[i][j]);
		}
		x[i] /= R[i][i];
	}
}

// ����������� ����� ����� QR-����������
template<typename variableType>
void Matrix<variableType>::QREigen(vector<variableType>& x) {
	Matrix<variableType> B(rows, cols);
	B = *this;
	Matrix<variableType> Q(rows, cols);
	Matrix<variableType> R(rows, cols);
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

// �������� ������� ����� LU-����������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::LUInverse() {
	// ���������� ��������� ������� 
	Matrix<variableType> E(rows, cols);

	// ������� ��� �����������
	Matrix<variableType> X(rows, cols);

	// � ����� �������� ������� ���� ����� LU-���������� �� �������� ������ X � E
	for (int i = 0; i < cols; i++)
		this->LUSolution(X[i], E[i]);

	// ������������� ������� X, ����� �������� ���������� ��� �������
	return X.transpose();
}

// �������� ������� ����� QR-����������
template<typename variableType>
Matrix<variableType> Matrix<variableType>::QRInverse() {
	// ���������� ������� ��� ����������
	Matrix<variableType> Q(rows, cols);
	Matrix<variableType> R(rows, cols);

	// �������� ������� Q � R
	this->QRDecomp(Q, R);

	// ������� �����������
	Matrix<variableType> X(rows, cols);

	// ������� ��� �������
	vector<variableType> y(rows, 0);
	vector<variableType> x(rows, 0);

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
template<typename variableType>
variableType Matrix<variableType>::norm() {
	variableType result = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++)
			result += pow((*this)[i][j], 2);
	}
	return sqrt(result);
}

// �������� ������� ������� ������ (Newton-Schulz-Hotelling)
template<typename variableType>
Matrix<variableType> Matrix<variableType>::ShultzInverse(double epsilon) {

	Matrix<variableType> A = *this;
	Matrix<variableType> E(rows, cols);
	Matrix<variableType> phi(rows, cols);
	Matrix<variableType> At = A.transpose();
	// �������� ��������� ������� X
	Matrix<variableType> U = At * (1.0 / (A * At).norm()); // ����� ���� ������� ��������

	for (int i = 0;; i++) {
		phi = E - A * U;
		if (phi.norm() < epsilon) break;
		U = U * (E + phi);
	}

	return U;
}

// ������� ���� ����� ����� ������-������
template<typename variableType>
void Matrix<variableType>::GaussSeidelSolution(vector<variableType>& x, const vector<variableType>& b, double epsilon) {
	Matrix<variableType> A = *this;
	vector<variableType> prev(x.size(), 0);

	// ��������� ������� �� ������������ ������������
	if (A.diagonal())
	{
		for (int k = 0; k < maxIteration; k++) {
			prev = x;
			for (int i = 0; i < rows; i++) {
				variableType var = 0;
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
template<typename variableType>
bool Matrix<variableType>::diagonal() {
	Matrix<variableType> A = *this;
	variableType sum;
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
template<typename variableType>
variableType Matrix<variableType>::RayleighEigen(double epsilon) {
	Matrix<variableType> A = *this;
	Matrix<variableType> E(rows, cols);
	vector<variableType> x(rows, 1);
	vector<variableType> y(rows, 0);
	variableType sum{};
	variableType lambda = scalarOp(A * x, x) / scalarOp(x, x);
	variableType p{};
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








/*
����� ���� ����� ��� ����������� ������ VectorT
*/






// operator[] ��� VectorT
template<typename variableType>
variableType& VectorT<variableType>::operator[](int i) {
	return v[i];
}

// operator[] ��� ������������ VectorT
template<typename variableType>
const variableType& VectorT<variableType>::operator[](int i) const {
	return v[i];
}

// ������ ������� VectorT
template<typename variableType>
int VectorT<variableType>::size() {
	return this->v.size();
}

// ��������� ������������������ ������� �� ������� - �����
template<typename variableType>
variableType VectorT<variableType>::operator* (vector<variableType>& v1) {
	variableType result(0);
	for (int i = 0; i < v.size(); i++)
		result += v[i] * v1[i];
	return result;
}

// ����� ����������������� ��������
template<typename variableType>
VectorT<variableType> VectorT<variableType>::operator+ (VectorT<variableType>& v1) {
	VectorT<variableType> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] + v1[i];
	return result;
}

// �������� ����������������� ��������
template<typename variableType>
VectorT<variableType> VectorT<variableType>::operator- (VectorT<variableType>& v1) {
	VectorT<variableType> result(v1.size());
	for (int i = 0; i < v1.size(); i++)
		result[i] = v[i] - v1[i];
	return result;
}

// �������� � ������ ������
template<typename variableType>
VectorT<variableType>  VectorT<variableType>::operator+ (variableType n) {
	VectorT<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] + n;
	}
	return result;
}

// �������� � ������ �����
template<typename variableType>
VectorT<variableType> operator+ (variableType n, VectorT<variableType>& v) {
	return v + n;
}

// ��������� ����� (�������� ������ ������)
template<typename variableType>
VectorT<variableType> VectorT<variableType>::operator- (variableType n) {
	VectorT<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] - n;
	}
	return result;
}

// ��������� �� ����� ������
template<typename variableType>
VectorT<variableType> VectorT<variableType>::operator* (variableType n) {
	VectorT<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n * v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n * v[i];
	}
	return result;
}

// ��������� �� ����� �����
template<typename variableType>
VectorT<variableType> operator* (variableType n, VectorT<variableType>& v) {
	return v * n;
}

// ������� �� ����� (�������� ������ ������)
template<typename variableType>
VectorT<variableType> VectorT<variableType>::operator/ (variableType n) {
	VectorT<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (abs(n / v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n / v[i];
	}
	return result;
}

// ����� ������������������ �������
template<typename variableType>
ostream& operator<< (ostream& out , VectorT<variableType>& vT) {
	for (double el : vT.v) {
		out << el << ' ';
	}
	return out;
}

// ���������� ������������������ �������
template<typename variableType>
istream& operator>> (istream& in, VectorT<variableType>& vT) {
	for (int i = 0; i < vT.size(); i++) {
		cout << " Input [" << i << "]" << " element : ";
		in >> vT[i];
	}
	return in;
}

// ��������� ��������� ��������
template<typename variableType>
Matrix<variableType> operator* (const vector<variableType>& v1, const VectorT<variableType>& v2) {
	Matrix<variableType> res(v1.size(), v1.size());
	for (int i = 0; i < v1.size(); i++) {
		for (int j = 0; j < v1.size(); j++)
			res[i][j] = v1[i] * v2[j];
	}
	return res;
}

// ���������������� �������� ������� - �������� ������ ������ VectorT
template<typename variableType>
VectorT<variableType> transpose(vector<variableType>& v1) {
	VectorT<variableType> v2(v1);
	return v2;
}









/*
����� ���� ����� ��� ����������� ��� ���������� �������� std::vector
*/









// ���������� ����� ������� ����� ��� ����������
template<typename variableType>
variableType vectorLength(const vector<variableType>& v) {
	variableType length = 0;
	for (int i = 0; i < v.size(); i++)
		length += std::pow(v[i], 2);
	return std::sqrt(length);
}

// ��������� ������� �� ����� �����
template<typename variableType>
vector<variableType> operator* (variableType n, const vector<variableType>& v) {
	vector<variableType> result(v.size(), 0);
	for (int i = 0; i < v.size(); i++) {
		if (abs(n * v[i]) < precision)
			result[i] = 0;
		else
			result[i] = n * v[i];
	}
	return result;
}

// ��������� ������� �� ����� ������
template<typename variableType>
vector<variableType> operator* (const vector<variableType>& v, variableType n) {
	return n * v;
}

// ��������� ����� �� �������
template<typename variableType>
vector<variableType> operator- (const vector<variableType>& v, variableType n) {
	vector<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] - n;
	}
	return result;
}


// �������� ������� � ������ �����
template<typename variableType>
vector<variableType> operator+ (variableType n, const vector<variableType>& v) {
	vector<variableType> result(v.size());
	for (int i = 0; i < v.size(); i++) {
		result[i] = v[i] + n;
	}
	return result;
}

// �������� ������� � ������ ������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>& v, variableType n) {
	return n + v;
}

// �������� ��������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>& v1, const vector<variableType>& v2) {
	vector<variableType> result(v1.size(), 0);
	for (int i = 0; i < v1.size(); i++)
		result[i] = v1[i] + v2[i];
	return result;
}

// ��������� ��������
template<typename variableType>
vector<variableType> operator- (const vector<variableType>& v1, const vector<variableType>& v2) {
	vector<variableType> result(v1.size(), 0);
	for (int i = 0; i < v1.size(); i++)
		result[i] = v1[i] - v2[i];
	return result;
}





// �������


template class VectorT<float>;
template class VectorT<double>;
template class VectorT<long double>;

template class Matrix<float>;
template class Matrix<double>;
template class Matrix<long double>;

template ostream& operator << (ostream& out, Matrix<float>& a);
template ostream& operator << (ostream& out, Matrix<double>& a);
template ostream& operator << (ostream& out, Matrix<long double>& a);

template ostream& operator << (ostream& out, VectorT<float>& a);
template ostream& operator << (ostream& out, VectorT<double>& a);
template ostream& operator << (ostream& out, VectorT<long double>& a);


template Matrix<float> operator* (float, Matrix<float>&);
template Matrix<double> operator* (double, Matrix<double>&);
template Matrix<long double> operator* (long double, Matrix<long double>&);


template vector<float> operator* (const vector<float>&, Matrix<float>&);
template vector<double> operator* (const vector<double>&, Matrix<double>&);
template vector<long double> operator* (const vector<long double>&, Matrix<long double>&);


template Matrix<float> operator+ (float, Matrix<float>&);
template Matrix<double> operator+ (double, Matrix<double>&);
template Matrix<long double> operator+ (long double, Matrix<long double>&);

template VectorT<float> operator+ (float, VectorT<float>&);
template VectorT<double> operator+ (double, VectorT<double>&);
template VectorT<long double> operator+ (long double, VectorT<long double>&);


template VectorT<float> operator* (float, VectorT<float>&);
template VectorT<double> operator* (double, VectorT<double>&);
template VectorT<long double> operator* (long double, VectorT<long double>&);

// ��� ���� ������ �������� �������� �� ��������
template Matrix<float> operator* (const vector<float>& v1, const VectorT<float>& v2);
template Matrix<double> operator* (const vector<double>& v1, const VectorT<double>& v2);
template Matrix<long double> operator* (const vector<long double>& v1, const VectorT<long double>& v2);

template VectorT<float> transpose(vector<float>& v1);
template VectorT<double> transpose(vector<double>& v1);
template VectorT<long double> transpose(vector<long double>& v1);