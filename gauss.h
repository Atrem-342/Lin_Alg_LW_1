#ifndef GAUSS_H
#define GAUSS_H

#include <stddef.h>

typedef enum GaussMode {
    GAUSS_MODE_NO_PIVOT = 0,
    GAUSS_MODE_PARTIAL_PIVOT = 1
} GaussMode;

typedef enum GaussStatus {
    GAUSS_OK = 0,
    GAUSS_INVALID_ARGUMENT,
    GAUSS_SINGULAR_MATRIX,
    GAUSS_ALLOCATION_FAILED
} GaussStatus;

GaussStatus gauss_solve(size_t n, const double *a, const double *b, double *x, GaussMode mode);

#endif
