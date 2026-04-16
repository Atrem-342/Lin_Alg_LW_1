#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct Matrix {
    size_t size;
    double *data;
} Matrix;

Matrix *matrix_create(size_t size);
void matrix_destroy(Matrix *matrix);
Matrix *matrix_copy(const Matrix *src);

int matrix_set(Matrix *matrix, size_t row, size_t col, double value);
int matrix_get(const Matrix *matrix, size_t row, size_t col, double *value);

void matrix_input(Matrix *matrix);
void matrix_print(const Matrix *matrix);

#endif
