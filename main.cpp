#include <iostream>
#include <iomanip>

// Наши библиотеки
#include "Matrix.h"
#include "LinearSystem.h"
#include "Eigen.h"
#include "Decomposition.h"
#include "Optimization.h"

using namespace std;

int main() {

	Matrix<double> E(4, 4);

	cin >> E;

	cout << E << '\n';

	Matrix<double> Q(4, 4);
	Matrix<double> R(4, 4);

	LUDecomp(E, Q, R);

	cout << Q << '\n';

	cout << R << '\n';

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