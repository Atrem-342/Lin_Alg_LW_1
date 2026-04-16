#include "analysis.h"

#include <math.h>

double vector_norm(size_t n, const double *vector) {
    double sum;
    size_t i;

    if (n == 0 || vector == NULL) {
        return 0.0;
    }

    sum = 0.0;
    for (i = 0; i < n; ++i) {
        sum += vector[i] * vector[i];
    }

    return sqrt(sum);
}

double residual_norm(const Matrix *matrix, const double *x, const double *b) {
    double diff_sum;
    double b_sum;
    size_t row;
    size_t col;

    if (matrix == NULL || matrix->data == NULL || x == NULL || b == NULL || matrix->size == 0) {
        return 0.0;
    }

    diff_sum = 0.0;
    b_sum = 0.0;

    for (row = 0; row < matrix->size; ++row) {
        double ax;
        double diff;

        ax = 0.0;
        for (col = 0; col < matrix->size; ++col) {
            ax += matrix->data[row * matrix->size + col] * x[col];
        }

        diff = ax - b[row];
        diff_sum += diff * diff;
        b_sum += b[row] * b[row];
    }

    if (b_sum == 0.0) {
        return sqrt(diff_sum);
    }

    return sqrt(diff_sum) / sqrt(b_sum);
}

double relative_error(size_t n, const double *exact, const double *approx) {
    double diff_sum;
    double exact_sum;
    size_t i;

    if (n == 0 || exact == NULL || approx == NULL) {
        return 0.0;
    }

    diff_sum = 0.0;
    exact_sum = 0.0;

    for (i = 0; i < n; ++i) {
        double diff;

        diff = approx[i] - exact[i];
        diff_sum += diff * diff;
        exact_sum += exact[i] * exact[i];
    }

    if (exact_sum == 0.0) {
        return sqrt(diff_sum);
    }

    return sqrt(diff_sum) / sqrt(exact_sum);
}
