#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <vector>
#include <iostream>
#include "Matrix.h"
#include <math.h>

using namespace std;

double f3_sphere(const vector<double>&);
double f3_hyper_ellipsoid(const vector<double>&);
double f3_test(const vector<double>&);
double f_diff(const vector<double>&, double (*)(const vector<double>&), int);

Matrix hessian(const vector<double>&, double (*)(const vector<double>&), vector<double>&);
vector<double> newton(const vector<double>&, double (*)(const vector<double>&));

vector<double> gradient(const vector<double>&, double (*)(const vector<double>&));
vector<double> gradientDescent(const vector<double>&, double (*)(const vector<double>&));

vector<double> BFGS(const vector<double>&, double (*)(const vector<double>&));


#endif