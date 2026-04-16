#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>

static int matrix_is_valid_index(const Matrix *matrix, size_t row, size_t col) {
    return matrix != NULL &&
           matrix->data != NULL &&
           row < matrix->size &&
           col < matrix->size;
}

Matrix *matrix_create(size_t size) {
    Matrix *matrix;

    if (size == 0) {
        return NULL;
    }

    matrix = (Matrix *)malloc(sizeof(Matrix));
    if (matrix == NULL) {
        return NULL;
    }

    matrix->size = size;
    matrix->data = (double *)calloc(size * size, sizeof(double));

    if (matrix->data == NULL) {
        free(matrix);
        return NULL;
    }

    return matrix;
}

void matrix_destroy(Matrix *matrix) {
    if (matrix == NULL) {
        return;
    }

    free(matrix->data);
    free(matrix);
}

Matrix *matrix_copy(const Matrix *src) {
    Matrix *copy;
    size_t i;
    size_t count;

    if (src == NULL || src->data == NULL || src->size == 0) {
        return NULL;
    }

    copy = matrix_create(src->size);
    if (copy == NULL) {
        return NULL;
    }

    count = src->size * src->size;
    for (i = 0; i < count; ++i) {
        copy->data[i] = src->data[i];
    }

    return copy;
}

int matrix_set(Matrix *matrix, size_t row, size_t col, double value) {
    if (!matrix_is_valid_index(matrix, row, col)) {
        return 0;
    }

    matrix->data[row * matrix->size + col] = value;
    return 1;
}

int matrix_get(const Matrix *matrix, size_t row, size_t col, double *value) {
    if (!matrix_is_valid_index(matrix, row, col) || value == NULL) {
        return 0;
    }

    *value = matrix->data[row * matrix->size + col];
    return 1;
}

void matrix_input(Matrix *matrix) {
    size_t row;
    size_t col;
    double value;

    if (matrix == NULL) {
        return;
    }

    for (row = 0; row < matrix->size; ++row) {
        for (col = 0; col < matrix->size; ++col) {
            printf("A[%zu][%zu] = ", row, col);
            if (scanf("%lf", &value) != 1) {
                return;
            }

            matrix_set(matrix, row, col, value);
        }
    }
}

void matrix_print(const Matrix *matrix) {
    size_t row;
    size_t col;
    double value;

    if (matrix == NULL || matrix->data == NULL) {
        printf("(пустая матрица)\n");
        return;
    }

    for (row = 0; row < matrix->size; ++row) {
        for (col = 0; col < matrix->size; ++col) {
            matrix_get(matrix, row, col, &value);
            printf("%10.6f ", value);
        }
        printf("\n");
    }
}
