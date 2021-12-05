#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include "Matrix.h"

// QR-разложение
template<typename variableType>
void QRDecomp(const Matrix<variableType>&, Matrix<variableType>&, Matrix<variableType>&);

// LU-разложение
template<typename variableType>
void LUDecomp(const Matrix<variableType>&, Matrix<variableType>&, Matrix<variableType>&);

// Разложение Холецкого
template<typename variableType>
void cholesky(const Matrix<variableType>&, Matrix<variableType>&);

#endif