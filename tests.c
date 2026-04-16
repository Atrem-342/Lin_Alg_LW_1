#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "analysis.h"
#include "gauss.h"
#include "generate_matrix.h"
#include "lu.h"
#include "matrix.h"

static void test_gauss_no_pivot(void) {
    const double a[] = {
        2.0, 1.0,
        1.0, 3.0
    };
    const double b[] = {
        5.0, 10.0
    };
    double x[2];

    assert(gauss_solve(2, a, b, x, GAUSS_MODE_NO_PIVOT) == GAUSS_OK);
    assert(fabs(x[0] - 1.0) < 1e-9);
    assert(fabs(x[1] - 3.0) < 1e-9);
}

static void test_gauss_partial_pivot(void) {
    const double a[] = {
        2.0, 1.0,
        1.0, 3.0
    };
    const double b[] = {
        5.0, 10.0
    };
    double x[2];

    assert(gauss_solve(2, a, b, x, GAUSS_MODE_PARTIAL_PIVOT) == GAUSS_OK);
    assert(fabs(x[0] - 1.0) < 1e-9);
    assert(fabs(x[1] - 3.0) < 1e-9);
}

static void test_lu_solve(void) {
    const double a[] = {
        2.0, 1.0,
        1.0, 3.0
    };
    const double b[] = {
        5.0, 10.0
    };
    double x[2];

    assert(lu_solve(2, a, b, x) == LU_OK);
    assert(fabs(x[0] - 1.0) < 1e-9);
    assert(fabs(x[1] - 3.0) < 1e-9);
}

static void test_lu_factored(void) {
    const double a[] = {
        2.0, 1.0,
        1.0, 3.0
    };
    const double b[] = {
        5.0, 10.0
    };
    double l[4];
    double u[4];
    double x[2];

    assert(lu_decompose(2, a, l, u) == LU_OK);
    assert(lu_solve_factored(2, l, u, b, x) == LU_OK);
    assert(fabs(x[0] - 1.0) < 1e-9);
    assert(fabs(x[1] - 3.0) < 1e-9);
}

static void test_hilbert_generation(void) {
    Matrix *matrix;

    matrix = generate_hilbert_matrix(3);
    assert(matrix != NULL);
    assert(fabs(matrix->data[0] - 1.0) < 1e-12);
    assert(fabs(matrix->data[1] - 0.5) < 1e-12);
    assert(fabs(matrix->data[4] - (1.0 / 3.0)) < 1e-12);
    matrix_destroy(matrix);
}

static void test_metrics(void) {
    Matrix *matrix;
    const double exact[] = {
        1.0, 3.0
    };
    const double b[] = {
        5.0, 10.0
    };

    matrix = matrix_create(2);
    assert(matrix != NULL);
    assert(matrix_set(matrix, 0, 0, 2.0));
    assert(matrix_set(matrix, 0, 1, 1.0));
    assert(matrix_set(matrix, 1, 0, 1.0));
    assert(matrix_set(matrix, 1, 1, 3.0));

    assert(residual_norm(matrix, exact, b) < 1e-12);
    assert(relative_error(2, exact, exact) < 1e-12);

    matrix_destroy(matrix);
}

int main(void) {
    test_gauss_no_pivot();
    test_gauss_partial_pivot();
    test_lu_solve();
    test_lu_factored();
    test_hilbert_generation();
    test_metrics();

    printf("Все тесты пройдены.\n");
    return 0;
}
