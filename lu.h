#ifndef LU_H
#define LU_H

#include <stddef.h>

typedef enum LUStatus {
    LU_OK = 0,
    LU_INVALID_ARGUMENT,
    LU_SINGULAR_MATRIX,
    LU_ALLOCATION_FAILED
} LUStatus;

LUStatus lu_decompose(size_t n, const double *a, double *l, double *u);
LUStatus forward_substitution(size_t n, const double *l, const double *b, double *y);
LUStatus backward_substitution(size_t n, const double *u, const double *y, double *x);
LUStatus lu_solve_factored(size_t n, const double *l, const double *u, const double *b, double *x);
LUStatus lu_solve(size_t n, const double *a, const double *b, double *x);

#endif
