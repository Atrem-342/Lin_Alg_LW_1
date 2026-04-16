#include "generate_matrix.h"

#include <stdlib.h>

static double random_value(double min_value, double max_value) {
    double t;

    t = (double)rand() / (double)RAND_MAX;
    return min_value + (max_value - min_value) * t;
}

void generate_matrix_seed(unsigned int seed) {
    srand(seed);
}

Matrix *generate_random_matrix(size_t size, double min_value, double max_value) {
    Matrix *matrix;
    size_t row;
    size_t col;
    double value;

    if (size == 0 || min_value > max_value) {
        return NULL;
    }

    matrix = matrix_create(size);
    if (matrix == NULL) {
        return NULL;
    }

    for (row = 0; row < size; ++row) {
        for (col = 0; col < size; ++col) {
            value = random_value(min_value, max_value);
            if (!matrix_set(matrix, row, col, value)) {
                matrix_destroy(matrix);
                return NULL;
            }
        }
    }

    return matrix;
}

double *generate_random_vector(size_t size, double min_value, double max_value) {
    double *vector;
    size_t i;

    if (size == 0 || min_value > max_value) {
        return NULL;
    }

    vector = (double *)malloc(size * sizeof(double));
    if (vector == NULL) {
        return NULL;
    }

    for (i = 0; i < size; ++i) {
        vector[i] = random_value(min_value, max_value);
    }

    return vector;
}

double *generate_constant_vector(size_t size, double value) {
    double *vector;
    size_t i;

    if (size == 0) {
        return NULL;
    }

    vector = (double *)malloc(size * sizeof(double));
    if (vector == NULL) {
        return NULL;
    }

    for (i = 0; i < size; ++i) {
        vector[i] = value;
    }

    return vector;
}

Matrix *generate_hilbert_matrix(size_t size) {
    Matrix *matrix;
    size_t row;
    size_t col;
    double value;

    if (size == 0) {
        return NULL;
    }

    matrix = matrix_create(size);
    if (matrix == NULL) {
        return NULL;
    }

    for (row = 0; row < size; ++row) {
        for (col = 0; col < size; ++col) {
            value = 1.0 / (double)(row + col + 1);
            if (!matrix_set(matrix, row, col, value)) {
                matrix_destroy(matrix);
                return NULL;
            }
        }
    }

    return matrix;
}

double *matrix_vector_multiply(const Matrix *matrix, const double *vector) {
    double *result;
    size_t row;
    size_t col;

    if (matrix == NULL || matrix->data == NULL || vector == NULL) {
        return NULL;
    }

    result = (double *)calloc(matrix->size, sizeof(double));
    if (result == NULL) {
        return NULL;
    }

    for (row = 0; row < matrix->size; ++row) {
        for (col = 0; col < matrix->size; ++col) {
            result[row] += matrix->data[row * matrix->size + col] * vector[col];
        }
    }

    return result;
}
