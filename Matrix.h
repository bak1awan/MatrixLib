#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;


// ��������������� ���������� ������ VectorT
template <class variableType>
class VectorT;








// ����� ������� ���������� ������ ������ Matrix









template <class variableType>
class Matrix {
public:
	// �����������, �� ��������� ���������� ��������� �������
	Matrix(int nRow = 1, int nCol = 1) : rows(nRow), cols(nCol) {
		arr = vector<vector<variableType>>(nRow, vector <variableType>(nCol, 0));
		for (int i = 0; i < nRow; i++) {
			arr[i][i] = 1;
		}
	};

	// ���������������� �������
	Matrix<variableType> transpose();

	// ������� �����
	Matrix<variableType> operator- ();

	// ���������� ��������� + � ��������
	Matrix<variableType> operator+ (const Matrix<variableType>&);

	// ���������� ��������� - � ��������
	Matrix<variableType> operator- (const Matrix<variableType>&);

	// ���������� ��������� * � ��������
	Matrix<variableType> operator* (const Matrix<variableType>&);

	// ��������� ������� �� ����� ������
	Matrix<variableType> operator* (variableType);

	// ��������� ������� �� ����� �����
	template <class T> friend Matrix<T> operator* (T, Matrix<T>&);

	// ��������� ������� �� ������ ������
	vector<variableType> operator* (const vector<variableType>&);

	// ��������� ������� �� ������ �����
	template <class T> friend vector<T> operator* (const vector<T>&, Matrix<T>&);

	// �������� ������� � ������ ������
	Matrix<variableType> operator+ (variableType);

	// �������� ������� � ������ �����
	template <class T> friend Matrix<T> operator+ (T, Matrix<T>&);

	// ���������� << ��� ������
	template <class T> friend ostream& operator << (std::ostream&, Matrix<T>&);

	// ���������� >> ��� ������
	template <class T> friend istream& operator >> (std::istream&, Matrix<T>&);

	// QR-����������
	void QRDecomp(Matrix<variableType>&, Matrix<variableType>&);

	// LU-����������
	void LUDecomp(Matrix<variableType>&, Matrix<variableType>&);

	// ������������ ����� QR-����������
	variableType LUDeterminant();

	// �������� ������� ����� LU-����������
	Matrix<variableType> LUInverse();

	// �������� ������� ����� QR-����������
	Matrix<variableType> QRInverse();

	// ������� ���� ������� ������ � �������
	void GaussSeidelSolution(vector<variableType>&, const vector<variableType>&, double = 0.001);

	// ����� �������
	variableType norm();

	// �������� ������� �� ������������ ������������
	bool diagonal();

	// �������� ������� ������� ������ (Newton-Schulz-Hotelling)
	Matrix<variableType> ShultzInverse(double = 0.001);

	// ������� ���� ����� LU-����������
	void LUSolution(vector<variableType>&, const vector<variableType>&);

	// ������� ���� ����� QR-����������
	void QRSolution(vector<variableType>&, const vector<variableType>&);

	void cholesky(Matrix<variableType>&);

	// ����������� ����� ����� QR-����������
	void QREigen(vector<variableType>&);

	// ����������� ����� ����� ����������� �����
	variableType RayleighEigen(double = 0.01);

	// ���������� ��������� [] ��� ������
	vector<variableType>& operator[] (const int i);

	// ���������� ��������� [] ��� ����������� ������
	const vector<variableType>& operator[] (const int i) const;
private:
	// ���� ������� �������������� ������ �������� ��������
	vector<vector<variableType>> arr;

	// ���������� �����
	int rows;

	// ���������� �������
	int cols;
};









// ����� ������� ���������� ��� ������ ������ ����������������� �������� VectorT









template <class variableType>
class VectorT {
public:
	// �����������
	VectorT(int size = 1) : v(size) {}

	// ����������� �����������
	VectorT(vector<variableType> v1) : v(v1) {}

	// ������ ������� VectorT
	int size();

	// ��������� ������������������ ������� �� ������� - �����
	variableType operator* (vector<variableType>&);

	// ����� ����������������� ��������
	VectorT<variableType> operator+ (VectorT<variableType>&);

	// �������� ����������������� ��������
	VectorT<variableType> operator- (VectorT<variableType>&);

	// �������� []
	variableType& operator[](int);
	const variableType& operator[](int i) const;

	// �������� � ������ ������
	VectorT<variableType> operator+ (variableType);

	// �������� � ������ �����
	template <class T> friend VectorT<T> operator+ (T, VectorT<T>&);

	// ��������� �����
	VectorT<variableType> operator- (variableType);

	// ��������� �� ����� ������
	VectorT<variableType> operator* (variableType);

	// ��������� �� ����� �����
	template <class T> friend VectorT<T> operator* (T, VectorT<T>&);

	// ������� �� ����� (�������� ������ ������)
	VectorT<variableType> operator/ (variableType);

	// ����� ������������������ �������
	template <class T> friend ostream& operator<< (ostream&, VectorT<T>&);

	// ���������� ������������������ �������
	template <class T> friend istream& operator>> (istream&, VectorT<T>&);

private:
	// ����, ���������� ������
	vector<variableType> v;
};








// ����� ������� ��������� ��� ���������� �������� std::vector








// ��������� ������� �� ����� �����
template<typename variableType>
vector<variableType> operator* (variableType, const vector<variableType>&);

// ��������� ������� �� ����� ������
template<typename variableType>
vector<variableType> operator* (const vector<variableType>&, variableType);

// �������� ������� � ������ �����
template<typename variableType>
vector<double> operator+ (variableType, const vector<variableType>&);

// �������� ������� � ������ ������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, variableType);

// ��������� ����� �� ������� (������ ������ ��������)
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, variableType);

// �������� ��������
template<typename variableType>
vector<variableType> operator+ (const vector<variableType>&, const vector<variableType>&);

// ��������� ��������
template<typename variableType>
vector<variableType> operator- (const vector<variableType>&, const vector<variableType>&);

// ���������� ����� �������
template<typename variableType>
variableType vectorLength(const vector<variableType>&);

// ��������� ������������ ��������
template<typename variableType>
variableType scalarOp(const vector<variableType>&, const vector<variableType>&);

// ���������������� �������� ������� - �������� ������ ������ VectorT
template<typename variableType>
VectorT<variableType> transpose(vector<variableType>&);

// ��������� ��������� ��������
template<typename variableType>
Matrix<variableType> operator* (const vector<variableType>&, const VectorT<variableType>&);
#endif