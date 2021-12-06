#include <iostream>
#include <iomanip>

// Наши библиотеки
#include "Matrix.h"
#include "LinearSystem.h"
#include "Eigen.h"
#include "Optimization.h"

using namespace std;

int main() {

	Matrix<double> E{ {2, 5, 7}, {6, 3, 4}, {5, -2, 3} };

	cout << E << '\n';

	Matrix<double> L = E.ShultzInverse();

	cout << setprecision(3) << L << '\n';

	// Matrix<double> E = v1 * v2;

	/*
	vector<vector<Matrix<double>>> v;
	v.resize(3);
	vector<Matrix<double>>::iterator it;
	for (int i = 0; i < v.size(); i++) {
		it = v[i].begin();
		for (int j = 0; j < v.size(); j++) {
			it = v[i].emplace(it, 3, 3);
		}
	}
	return 0;
	*/

}