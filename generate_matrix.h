#ifndef GENERATE_MATRIX_H
#define GENERATE_MATRIX_H

#include <stddef.h>
#include "matrix.h"

void generate_matrix_seed(unsigned int seed);

Matrix *generate_random_matrix(size_t size, double min_value, double max_value);
double *generate_random_vector(size_t size, double min_value, double max_value);
double *generate_constant_vector(size_t size, double value);
Matrix *generate_hilbert_matrix(size_t size);
double *matrix_vector_multiply(const Matrix *matrix, const double *vector);

#endif
