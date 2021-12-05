#ifndef EIGEN_H
#define EIGEN_H

#include "Matrix.h"

// Собственные числа через QR-разложение
template <typename variableType>
void QREigen(const Matrix<variableType>&, vector<variableType>&);

// Собственное число через соотношение Рэлея
template <typename variableType>
variableType RayleighEigen(const Matrix<variableType>&, variableType);

#endif