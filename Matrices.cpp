#include <iostream >
#include "Matrix.h"
#include <iomanip>
#include "Optimization.h"

using namespace std;

int main() {
	
	Matrix A(3, 3);
	
	cout << " Input matrix A: " << endl;
	cin >> A;

	Matrix L(3, 3);

	A.cholesky(L);

	cout << L;

	Matrix LT = L.transpose();

	cout << L * LT;

	/*

	Matrix B(A);

	Matrix C = A * B;

	cout << "A: \n";

	cout << A;

	cout << "B: \n";

	cout << B;

	cout << "C: \n";

	cout << C;

	*/
	
	


	// Максимальное собственное значение через метод Рэлея
	/*
	double lambda = A.RayleighEigen(0.0001);
	cout << lambda;
	*/
	


	// Решение СЛАУ через метод Гаусса-Зейделя
	/*
	vector<double> b{ 5, 4, 3, 2 };
	vector<double> x{ 1, 1, 1, 1 };
	A.GaussSeidelSolution(x, b, 0.001);
	for (int i = 0; i < x.size(); i++)
		cout << x[i] << '\t';
	*/

	// Обратная матрица через метод Шульца
	/*
	cout << setprecision(3) << A.ShultzInverse(0.0001);
	*/
	
	// Обратная матрица через LU-разложение
	/*
	Matrix X = A.LUInverse();
	cout << setprecision(3) << X;
	*/

	// Обратная матрица через QR-разложение
	/*
	Matrix X = A.QRInverse();
	cout << setprecision(3) << X;
	*/
	
	
	
	/*
	cout << " matrix A: " << endl;
	cout << A;
	*/

	/*
	cout << "A + B:" << endl;
	C = A + B;
	C.print();
	*/

	/*
	cout << "A - B:" << endl;
	C = A - B;
	C.print();
	*/

	/*
	cout << "A * B:" << endl;
	Matrix C = A * B;
	cout << C;
	*/

	/*
	int n = 4;
	cout << "A * " << n << ":\n";
	Matrix D = A * n;
	cout << D;
	*/

	/*
	cout << "\nTransposed A:\n";
	Matrix C = A.transpose();
	cout << C;
	*/

	// QR-разложение
	/*
	Matrix Q(A.rows, A.cols);

	cout << setprecision(3) << "Matrix Q:\n";
	cout << Q;

	Matrix R(A.rows, A.cols);

	cout << "Matrix R:\n";
	cout << setprecision(3) << R;

	A.QRDecomp(Q, R);

	cout << "Matrix Q after QRDecomp:\n";
	cout << setprecision(3) << Q;

	cout << "Matrix R after QRDecomp:\n";
	cout << setprecision(3) << R;

	cout << "Q * R:\n";
	cout << setprecision(3) << Q * R;
	*/
	


	// LU-разложение
	/*
	Matrix L(A.rows, A.cols);
	Matrix U(A.rows, A.cols);
	A.LUDecomp(L, U);

	cout << setprecision(3) << "Matrix A:\n";
	cout << setprecision(3) << A;

	cout << setprecision(3) << "Matrix L:\n";
	cout << setprecision(3) << L;

	cout << "Matrix U:\n";
	cout << setprecision(3) << U;
	*/
	
	

	// Определитель через LU-разложение
	/*
	double determinant = A.LUDeterminant();
	cout << determinant;
	*/

	// Решение СЛАУ через LU-разложение
	/*
	vector<double> b{ 0, 3, 1 };
	vector<double> x{ 0, 0, 0 };
	A.LUSolution(x, b);
	for (int i = 0; i < x.size(); i++)
		cout << setprecision(3) << x[i] << '\t';
	*/
	
	/*
	vector<double> x{ 0, 0 };
		A.QREigen(x);
		for (int i = 0; i < x.size(); i++)
			cout << x[i] << '\t';
	*/
	

	// Решение СЛАУ через QR-разложение
	/*
	vector<double> b{ 1, 2, 3 };
	vector<double> x(b.size(), 0);
	A.QRSolution(x, b);
	for (int i = 0; i < x.size(); i++)
		cout << setprecision(3) << x[i] << '\t';
	*/
	
	// Градиентный спуск
	
	vector<double> x1{ 4.9, 4.9 };
	x1 = gradientDescent(x1, f3_hyper_ellipsoid);
	for (int i = 0; i < x1.size(); i++)
		cout << x1[i] << '\t';

	cout << '\n';
	


	// Ньютон
	
	vector<double> x2{ 4.9, 4.9 };
	x2 = newton(x2, f3_hyper_ellipsoid);
	for (int i = 0; i < x2.size(); i++)
		cout << x2[i] << '\t';	
	
	cout << '\n';

	// BFGS
	
	vector<double> x3{ 4.9, 4.9 };
	x3 = BFGS(x3, f3_hyper_ellipsoid);
	for (int i = 0; i < x3.size(); i++)
		cout << x3[i] << '\t';
	return 0;
	

}