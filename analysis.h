#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stddef.h>

#include "matrix.h"

double vector_norm(size_t n, const double *vector);
double residual_norm(const Matrix *matrix, const double *x, const double *b);
double relative_error(size_t n, const double *exact, const double *approx);

#endif
