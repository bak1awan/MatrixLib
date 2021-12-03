#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

using namespace std;

class VectorT;

class Matrix {
public:
	// �����������
	int rows;
	int cols;
	Matrix(int nRow = 1, int nCol = 1);

	// ���������������� �������
	Matrix transpose();

	// ������� �����
	Matrix operator- ();

	// ���������� ��������� +
	Matrix operator+ (const Matrix&);

	// ���������� ��������� -
	Matrix operator- (const Matrix&);

	// ���������� ��������� * � ��������
	Matrix operator* (const Matrix&);

	// ���������� ��������� * � ������ ������
	friend Matrix operator* (Matrix&, double);

	// ���������� ��������� * � ������ �����
	friend Matrix operator* (double, Matrix&);

	// ���������� ��������� * � �������� ������
	friend vector<double> operator* (Matrix&, const vector<double>&);

	// ���������� ��������� * � �������� �����
	friend vector<double> operator* (const vector<double>&, Matrix&);

	// ���������� ��������� + � ������ ������
	friend Matrix operator+ (Matrix&, double);

	// ���������� ��������� + � ������ �����
	friend Matrix operator+ (double, Matrix&);

	// ���������� <<
	friend std::ostream& operator << (std::ostream&, const Matrix&);

	// ���������� >>
	friend std::istream& operator >> (std::istream&, Matrix&);

	// QR-����������
	void QRDecomp(Matrix&, Matrix&);

	// LU-����������
	void LUDecomp(Matrix&, Matrix&);

	// ������������ ����� QR-����������
	double LUDeterminant();

	// �������� ������� ����� LU-����������
	Matrix LUInverse();

	// �������� ������� ����� QR-����������
	Matrix QRInverse();

	// ������� ���� ������� ������ � �������
	void GaussSeidelSolution(vector<double>&, const vector<double>&, double = 0.001);

	// ����� �������
	double norm();

	// �������� ������� �� ������������ ������������
	bool diagonal();

	// �������� ������� ������� ������ (Newton-Schulz-Hotelling)
	Matrix ShultzInverse(double = 0.001);

	// ������� ���� ����� LU-����������
	void LUSolution(vector<double>&, const vector<double>&);

	// ������� ���� ����� QR-����������
	void QRSolution(vector<double>&, const vector<double>&);

	void cholesky(Matrix&);

	// ����������� ����� ����� QR-����������
	void QREigen(vector<double>&);

	// ����������� ����� ����� ����������� �����
	double RayleighEigen(double = 0.01);

	// ���������� ��������� [] ��� ������
	vector<double>& operator[] (const int i);

	const vector<double>& operator[] (const int i) const;

	// ���� ��� �������������
	vector<vector<double>> arr;
};

class VectorT {
public:
	VectorT(int size = 1) : v(size) {}

	// ����������� �����������
	VectorT(vector<double>& v1) : v(v1) {}

	// ����������� �����������
	VectorT(vector<double>&& v1) : v(v1) {}

	// ����������� ������������ ������������
	VectorT operator= (vector<double>& v1) {
		this->v = v1;
	}

	// ����������� ������������ ������������
	VectorT operator= (vector<double>&& v1) {
		this->v = v1;
	}

	int size();

	// ��������� ������������������ ������� �� ������� - �����
	double operator* (vector<double>&);

	// ����� ����������������� ��������
	VectorT operator+ (VectorT&);

	// �������� ����������������� ��������
	VectorT operator- (VectorT&);

	// �������� []
	double& operator[](int);
	const double& operator[](int i) const;

	// �������� � ������
	VectorT operator+ (double);
	friend VectorT operator+ (double, VectorT&);

	// ��������� �����
	VectorT operator- (double);

	// ��������� �� �����
	VectorT operator* (double);
	friend VectorT operator* (double, VectorT&);

	// ������� �� �����
	VectorT operator/ (double);

	// ����� ������������������ �������
	friend ostream& operator<< (ostream&, VectorT&);

	// ���������� ������������������ �������
	friend istream& operator>> (istream&, VectorT&);

private:
	// ����, ���������� ������
	vector<double> v;
};

// ��������� ������� �� �����
vector<double> operator* (double, const vector<double>&);
vector<double> operator* (const vector<double>&, double);

// �������� ������� � ������
vector<double> operator+ (double, const vector<double>&);
vector<double> operator+ (const vector<double>&, double);

vector<double> operator- (const vector<double>&, double);

// �������� ��������
vector<double> operator+ (const vector<double>&, const vector<double>&);

// ��������� ��������
vector<double> operator- (const vector<double>&, const vector<double>&);

// ���������� ����� �������
double vectorLength(const vector<double>&);

// ��������� ������������ ��������
Matrix operator*(const vector<double>&, const vector<double>&);

// ��������� ������������ ��������
double scalarOp(const vector<double>&, const vector<double>&);

VectorT transpose(vector<double>&);

// ��������� ��������� �������
Matrix generateE(int);

// ��������� �������� ������� �� ����������������� - �������
Matrix operator* (const vector<double>&, const VectorT&);

#endif