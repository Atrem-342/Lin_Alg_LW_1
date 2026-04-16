#include "lu.h"

#include <math.h>
#include <stdlib.h>

static int is_near_zero(double value) {
    return fabs(value) < 1e-18;
}

static void set_identity(size_t n, double *m) {
    size_t row;
    size_t col;

    for (row = 0; row < n; ++row) {
        for (col = 0; col < n; ++col) {
            m[row * n + col] = (row == col) ? 1.0 : 0.0;
        }
    }
}

LUStatus lu_decompose(size_t n, const double *a, double *l, double *u) {
    size_t row;
    size_t col;
    size_t k;

    if (n == 0 || a == NULL || l == NULL || u == NULL) {
        return LU_INVALID_ARGUMENT;
    }

    set_identity(n, l);

    for (row = 0; row < n * n; ++row) {
        u[row] = 0.0;
    }

    for (k = 0; k < n; ++k) {
        for (col = k; col < n; ++col) {
            double sum = 0.0;

            for (row = 0; row < k; ++row) {
                sum += l[k * n + row] * u[row * n + col];
            }

            u[k * n + col] = a[k * n + col] - sum;
        }

        if (is_near_zero(u[k * n + k])) {
            return LU_SINGULAR_MATRIX;
        }

        for (row = k + 1; row < n; ++row) {
            double sum = 0.0;

            for (col = 0; col < k; ++col) {
                sum += l[row * n + col] * u[col * n + k];
            }

            l[row * n + k] = (a[row * n + k] - sum) / u[k * n + k];
        }
    }

    return LU_OK;
}

LUStatus forward_substitution(size_t n, const double *l, const double *b, double *y) {
    size_t row;
    size_t col;

    if (n == 0 || l == NULL || b == NULL || y == NULL) {
        return LU_INVALID_ARGUMENT;
    }

    for (row = 0; row < n; ++row) {
        double sum = b[row];

        for (col = 0; col < row; ++col) {
            sum -= l[row * n + col] * y[col];
        }

        if (is_near_zero(l[row * n + row])) {
            return LU_SINGULAR_MATRIX;
        }

        y[row] = sum / l[row * n + row];
    }

    return LU_OK;
}

LUStatus backward_substitution(size_t n, const double *u, const double *y, double *x) {
    size_t row;
    size_t col;

    if (n == 0 || u == NULL || y == NULL || x == NULL) {
        return LU_INVALID_ARGUMENT;
    }

    for (row = n; row-- > 0;) {
        double sum = y[row];

        for (col = row + 1; col < n; ++col) {
            sum -= u[row * n + col] * x[col];
        }

        if (is_near_zero(u[row * n + row])) {
            return LU_SINGULAR_MATRIX;
        }

        x[row] = sum / u[row * n + row];
    }

    return LU_OK;
}

LUStatus lu_solve_factored(size_t n, const double *l, const double *u, const double *b, double *x) {
    double *y;
    LUStatus status;

    if (n == 0 || l == NULL || u == NULL || b == NULL || x == NULL) {
        return LU_INVALID_ARGUMENT;
    }

    y = (double *)malloc(n * sizeof(double));
    if (y == NULL) {
        return LU_ALLOCATION_FAILED;
    }

    status = forward_substitution(n, l, b, y);
    if (status == LU_OK) {
        status = backward_substitution(n, u, y, x);
    }

    free(y);
    return status;
}

LUStatus lu_solve(size_t n, const double *a, const double *b, double *x) {
    double *l;
    double *u;
    LUStatus status;

    if (n == 0 || a == NULL || b == NULL || x == NULL) {
        return LU_INVALID_ARGUMENT;
    }

    l = (double *)malloc(n * n * sizeof(double));
    if (l == NULL) {
        return LU_ALLOCATION_FAILED;
    }

    u = (double *)malloc(n * n * sizeof(double));
    if (u == NULL) {
        free(l);
        return LU_ALLOCATION_FAILED;
    }

    status = lu_decompose(n, a, l, u);
    if (status == LU_OK) {
        status = lu_solve_factored(n, l, u, b, x);
    }

    free(l);
    free(u);
    return status;
}
