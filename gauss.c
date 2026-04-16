#include "gauss.h"
#include <math.h>
#include <stdlib.h>

static int is_near_zero(double value) {
    return fabs(value) < 1e-18;
}

// FOR BOTH

static double *matrix_copy(size_t n, const double *src) {
    double *dst;
    size_t count;

    if (src == NULL) {
        return NULL;
    }

    count = n * n;
    dst = (double *)malloc(count * sizeof(double));
    if (dst == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < count; ++i) {
        dst[i] = src[i];
    }

    return dst;
}

static double *vector_copy(size_t n, const double *src) {
    double *dst;

    if (src == NULL) {
        return NULL;
    }

    dst = (double *)malloc(n * sizeof(double));
    if (dst == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < n; ++i) {
        dst[i] = src[i];
    }

    return dst;
}

// FOR PARTIAL PIVOT

static void swap_rows(size_t n, double *a, size_t row1, size_t row2) {
    size_t col;

    if (row1 == row2){
        return;
    }

    for (col = 0; col < n; ++col){
        double tmp = a[row1 *n +col];
        a[row1 * n + col] = a[row2 * n + col];
        a[row2 * n + col] = tmp;
    }
}

static void swap_values(double *x, double *y) {
    double tmp = *x;
    *x = *y;
    *y = tmp;
}

// FOR NO PIVOT 

static GaussStatus forward_elimination(size_t n, double *a, double *rhs) {
    size_t pivot_row;
    size_t row;
    size_t col;

    for (pivot_row = 0; pivot_row + 1 < n; ++pivot_row) {
        double pivot = a[pivot_row * n + pivot_row];
        if (is_near_zero(pivot)) {
            return GAUSS_SINGULAR_MATRIX;
        }

        for (row = pivot_row + 1; row < n; ++row) {
            double factor = a[row * n + pivot_row] / pivot;
            a[row * n + pivot_row] = 0.0;

            for (col = pivot_row + 1; col < n; ++col) {
                a[row * n + col] -= factor * a[pivot_row * n + col];
            }

            rhs[row] -= factor * rhs[pivot_row];
        }
    }

    if (is_near_zero(a[(n - 1) * n + (n - 1)])) {
        return GAUSS_SINGULAR_MATRIX;
    }

    return GAUSS_OK;
}

// FOR PARTIAL PIVOT

static GaussStatus forward_elimination_partial_pivot(size_t n, double *a, double *rhs) {
    size_t pivot_col;
    size_t row;
    size_t col;

    for (pivot_col = 0; pivot_col + 1 < n; ++pivot_col) {
        size_t pivot_row = pivot_col;
        double max_value = fabs( a[pivot_col * n + pivot_col] );

        for (row = pivot_col +1; row < n; ++row) {
            double current = fabs( a[row * n + pivot_col]);
            if (current > max_value) {
                max_value = current;
                pivot_row = row;
            }
        }

        if (is_near_zero(max_value)) {
            return GAUSS_SINGULAR_MATRIX;
        }

        if (pivot_row != pivot_col) {
            swap_rows(n, a, pivot_col, pivot_row);
            swap_values( &rhs[pivot_col], &rhs[pivot_row]);
        }

        {
            double pivot = a[pivot_col * n + pivot_col];

            if (is_near_zero(pivot)) {
                return GAUSS_SINGULAR_MATRIX;
            }

            for (row = pivot_col + 1; row < n; ++row) {
                double factor = a[row * n + pivot_col] / pivot;
                a[row * n + pivot_col] = 0.0;

                for (col = pivot_col + 1; col < n; ++col){
                    a[row * n + col] -= factor* a[pivot_col * n + col];
                }

                rhs[row] -= factor * rhs[pivot_col];
            }
        }
    }

    if (is_near_zero( a[(n-1) * n + (n-1)]) ) {
        return GAUSS_SINGULAR_MATRIX;
    }
    
    return GAUSS_OK;
}

// FOR BOTH

static GaussStatus back_substitution(size_t n, const double *a, const double *rhs, double *x) {
    size_t row;

    for (row = n; row-- > 0;) {
        double sum = rhs[row];
        size_t col;

        for (col = row + 1; col < n; ++col) {
            sum -= a[row * n + col] * x[col];
        }

        if (is_near_zero(a[row * n + row])) {
            return GAUSS_SINGULAR_MATRIX;
        }

        x[row] = sum / a[row * n + row];
    }

    return GAUSS_OK;
}

GaussStatus gauss_solve(size_t n, const double *a, const double *b, double *x, GaussMode mode) {
    double *a_work;
    double *b_work;
    GaussStatus status;

    if (n == 0 || a == NULL || b == NULL || x == NULL) {
        return GAUSS_INVALID_ARGUMENT;
    }

    a_work = matrix_copy(n, a);
    if (a_work == NULL) {
        return GAUSS_ALLOCATION_FAILED;
    }

    b_work = vector_copy(n, b);
    if (b_work == NULL) {
        free(a_work);
        return GAUSS_ALLOCATION_FAILED;
    }
    switch (mode){
        case GAUSS_MODE_NO_PIVOT:
            status = forward_elimination(n, a_work, b_work);
            break;

        case GAUSS_MODE_PARTIAL_PIVOT:
            status = forward_elimination_partial_pivot(n, a_work, b_work);
            break;
        default:
            free(a_work);
            free(b_work);
            return GAUSS_INVALID_ARGUMENT;   
    }

    if (status == GAUSS_OK) {
        status = back_substitution(n, a_work, b_work, x);
    }

    free(a_work);
    free(b_work);
    return status;
} 
