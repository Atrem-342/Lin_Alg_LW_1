#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "analysis.h"
#include "clock_get_time.h"
#include "gauss.h"
#include "generate_matrix.h"
#include "lu.h"
#include "matrix.h"

typedef struct GaussTask {
    size_t n;
    const double *a;
    const double *b;
    double *x;
    GaussMode mode;
    GaussStatus status;
} GaussTask;

typedef struct LUTask {
    size_t n;
    const double *a;
    const double *b;
    double *x;
    LUStatus status;
} LUTask;

typedef struct LUDecomposeTask {
    size_t n;
    const double *a;
    double *l;
    double *u;
    LUStatus status;
} LUDecomposeTask;

typedef struct LUFactoredTask {
    size_t n;
    const double *l;
    const double *u;
    const double *b;
    double *x;
    LUStatus status;
} LUFactoredTask;

static void run_gauss_task(void *arg) {
    GaussTask *task = (GaussTask *)arg;

    task->status = gauss_solve(task->n, task->a, task->b, task->x, task->mode);
}

static void run_lu_task(void *arg) {
    LUTask *task = (LUTask *)arg;

    task->status = lu_solve(task->n, task->a, task->b, task->x);
}

static void run_lu_decompose_task(void *arg) {
    LUDecomposeTask *task = (LUDecomposeTask *)arg;

    task->status = lu_decompose(task->n, task->a, task->l, task->u);
}

static void run_lu_factored_task(void *arg) {
    LUFactoredTask *task = (LUFactoredTask *)arg;

    task->status = lu_solve_factored(task->n, task->l, task->u, task->b, task->x);
}

static double *allocate_vector(size_t n) {
    return (double *)calloc(n, sizeof(double));
}

static int read_matrix(Matrix *matrix) {
    size_t row;
    size_t col;
    double value;

    if (matrix == NULL || matrix->data == NULL) {
        return 0;
    }

    for (row = 0; row < matrix->size; ++row) {
        for (col = 0; col < matrix->size; ++col) {
            printf("A[%zu][%zu] = ", row, col);
            if (scanf("%lf", &value) != 1) {
                return 0;
            }

            if (!matrix_set(matrix, row, col, value)) {
                return 0;
            }
        }
    }

    return 1;
}

static int read_vector(size_t n, double *vector) {
    size_t i;

    if (vector == NULL) {
        return 0;
    }

    for (i = 0; i < n; ++i) {
        printf("b[%zu] = ", i);
        if (scanf("%lf", &vector[i]) != 1) {
            return 0;
        }
    }

    return 1;
}

static void print_vector(size_t n, const double *vector) {
    size_t i;

    for (i = 0; i < n; ++i) {
        printf("x[%zu] = %.10f\n", i, vector[i]);
    }
}

static void make_diagonally_dominant(Matrix *matrix) {
    size_t row;
    size_t col;

    if (matrix == NULL || matrix->data == NULL) {
        return;
    }

    for (row = 0; row < matrix->size; ++row) {
        double off_diagonal_sum = 0.0;

        for (col = 0; col < matrix->size; ++col) {
            if (row != col) {
                off_diagonal_sum += fabs(matrix->data[row * matrix->size + col]);
            }
        }

        matrix->data[row * matrix->size + row] += off_diagonal_sum + 1.0;
    }
}

static Matrix *create_random_system_matrix(size_t n) {
    Matrix *matrix;

    matrix = generate_random_matrix(n, -1.0, 1.0);
    if (matrix == NULL) {
        return NULL;
    }

    make_diagonally_dominant(matrix);
    return matrix;
}

static void print_main_menu(void) {
    printf("\nМеню:\n");
    printf("1 - Ручной ввод и решение одной системы\n");
    printf("2 - Эксперимент: случайные матрицы и время\n");
    printf("3 - Эксперимент: несколько правых частей k = 1, 10, 100\n");
    printf("4 - Эксперимент: матрица Гильберта и точность\n");
    printf("0 - Выход\n");
    printf("Ваш выбор: ");
}

static void print_method_menu(void) {
    printf("\nМетод:\n");
    printf("1 - Гаусс без выбора главного элемента\n");
    printf("2 - Гаусс с частичным выбором главного элемента\n");
    printf("3 - LU-разложение\n");
    printf("Ваш выбор: ");
}

static void print_rhs_table_header(void) {
    printf("%-5s %-5s %-22s %-14s %-14s %-14s %-14s %-6s\n",
           "n", "k", "method", "decompose_sec", "solve_sec", "total_sec", "avg_sec", "status");
}

static void print_rhs_table_row(size_t n,
                                size_t k,
                                const char *method,
                                double decompose_sec,
                                double solve_sec,
                                double total_sec,
                                double avg_sec,
                                int status) {
    printf("%-5zu %-5zu %-22s %-14.6e %-14.6e %-14.6e %-14.6e %-6d\n",
           n, k, method, decompose_sec, solve_sec, total_sec, avg_sec, status);
}

static void print_solution_block(const char *title,
                                 double elapsed,
                                 int status,
                                 const Matrix *matrix,
                                 const double *x,
                                 const double *b,
                                 const double *x_exact) {
    printf("%s\n", title);
    printf("Время: %.6f сек\n", elapsed);
    printf("Статус: %d\n", status);

    if (status == 0 && matrix != NULL && x != NULL && b != NULL) {
        print_vector(matrix->size, x);
        printf("Невязка: %.6e\n", residual_norm(matrix, x, b));
        if (x_exact != NULL) {
            printf("Относительная ошибка: %.6e\n", relative_error(matrix->size, x_exact, x));
        }
    }
}

static void run_manual_case(void) {
    size_t n;
    int method;
    Matrix *a;
    double *b;
    double *x;

    printf("Введите размер системы n: ");
    if (scanf("%zu", &n) != 1 || n == 0) {
        printf("Ошибка ввода размера.\n");
        return;
    }

    a = matrix_create(n);
    b = allocate_vector(n);
    x = allocate_vector(n);

    if (a == NULL || b == NULL || x == NULL) {
        printf("Ошибка выделения памяти.\n");
        matrix_destroy(a);
        free(b);
        free(x);
        return;
    }

    printf("\nВведите матрицу A:\n");
    if (!read_matrix(a)) {
        printf("Ошибка ввода матрицы.\n");
        matrix_destroy(a);
        free(b);
        free(x);
        return;
    }

    printf("\nВведите вектор b:\n");
    if (!read_vector(n, b)) {
        printf("Ошибка ввода вектора.\n");
        matrix_destroy(a);
        free(b);
        free(x);
        return;
    }

    print_method_menu();
    if (scanf("%d", &method) != 1) {
        printf("Ошибка ввода метода.\n");
        matrix_destroy(a);
        free(b);
        free(x);
        return;
    }

    if (method == 1 || method == 2) {
        GaussTask task;
        double elapsed;

        task.n = n;
        task.a = a->data;
        task.b = b;
        task.x = x;
        task.mode = (method == 1) ? GAUSS_MODE_NO_PIVOT : GAUSS_MODE_PARTIAL_PIVOT;
        task.status = GAUSS_OK;

        elapsed = measure_seconds(run_gauss_task, &task);
        print_solution_block("Решение методом Гаусса:", elapsed, task.status, a, x, b, NULL);
    } else if (method == 3) {
        LUTask task;
        double elapsed;

        task.n = n;
        task.a = a->data;
        task.b = b;
        task.x = x;
        task.status = LU_OK;

        elapsed = measure_seconds(run_lu_task, &task);
        print_solution_block("Решение методом LU:", elapsed, task.status, a, x, b, NULL);
    } else {
        printf("Неизвестный метод.\n");
    }

    matrix_destroy(a);
    free(b);
    free(x);
}

static void run_random_timing_experiment(void) {
    const size_t sizes[] = {100, 200, 500, 1000};
    const size_t count = sizeof(sizes) / sizeof(sizes[0]);
    size_t index;

    printf("n,method,status,time_sec,residual,relative_error\n");

    for (index = 0; index < count; ++index) {
        size_t n = sizes[index];
        Matrix *a;
        double *x_exact;
        double *b;
        double *x;

        a = create_random_system_matrix(n);
        x_exact = generate_constant_vector(n, 1.0);
        b = (a != NULL && x_exact != NULL) ? matrix_vector_multiply(a, x_exact) : NULL;
        x = allocate_vector(n);

        if (a == NULL || x_exact == NULL || b == NULL || x == NULL) {
            printf("%zu,error,%d,0.000000e+00,-1.000000e+00,-1.000000e+00\n", n, -1);
            matrix_destroy(a);
            free(x_exact);
            free(b);
            free(x);
            continue;
        }

        {
            GaussTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.mode = GAUSS_MODE_NO_PIVOT;
            task.status = GAUSS_OK;

            elapsed = measure_seconds(run_gauss_task, &task);
            residual = (task.status == GAUSS_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == GAUSS_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,gauss_no_pivot,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        {
            GaussTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.mode = GAUSS_MODE_PARTIAL_PIVOT;
            task.status = GAUSS_OK;

            elapsed = measure_seconds(run_gauss_task, &task);
            residual = (task.status == GAUSS_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == GAUSS_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,gauss_partial_pivot,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        {
            LUTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.status = LU_OK;

            elapsed = measure_seconds(run_lu_task, &task);
            residual = (task.status == LU_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == LU_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,lu,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        matrix_destroy(a);
        free(x_exact);
        free(b);
        free(x);
    }
}

static void run_multi_rhs_experiment(void) {
    const size_t n = 200;
    const size_t rhs_counts[] = {1, 10, 100};
    const size_t count = sizeof(rhs_counts) / sizeof(rhs_counts[0]);
    size_t index;

    printf("Используется размер n = %zu\n", n);
    print_rhs_table_header();

    for (index = 0; index < count; ++index) {
        size_t k = rhs_counts[index];
        Matrix *a;
        double *l;
        double *u;
        LUDecomposeTask decompose_task;
        double decompose_time;
        double gauss_total;
        double lu_total;
        int gauss_status;
        int lu_status;
        size_t rhs;

        a = create_random_system_matrix(n);
        l = (double *)malloc(n * n * sizeof(double));
        u = (double *)malloc(n * n * sizeof(double));

        if (a == NULL || l == NULL || u == NULL) {
            printf("%zu,%zu,error,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-1\n", n, k);
            matrix_destroy(a);
            free(l);
            free(u);
            continue;
        }

        decompose_task.n = n;
        decompose_task.a = a->data;
        decompose_task.l = l;
        decompose_task.u = u;
        decompose_task.status = LU_OK;

        decompose_time = measure_seconds(run_lu_decompose_task, &decompose_task);
        if (decompose_task.status != LU_OK) {
            printf("%zu,%zu,lu_decompose,%d,%.6e,0.000000e+00,0.000000e+00,%.6e\n",
                   n, k, decompose_task.status, decompose_time, decompose_time);
            matrix_destroy(a);
            free(l);
            free(u);
            continue;
        }

        gauss_total = 0.0;
        lu_total = 0.0;
        gauss_status = GAUSS_OK;
        lu_status = LU_OK;

        for (rhs = 0; rhs < k; ++rhs) {
            double scale;
            double *x_exact;
            double *b;
            double *x_gauss;
            double *x_lu;
            GaussTask gauss_task;
            LUFactoredTask lu_task;
            double elapsed;

            scale = 1.0 + 0.1 * (double)rhs;
            x_exact = generate_constant_vector(n, scale);
            b = (x_exact != NULL) ? matrix_vector_multiply(a, x_exact) : NULL;
            x_gauss = allocate_vector(n);
            x_lu = allocate_vector(n);

            if (x_exact == NULL || b == NULL || x_gauss == NULL || x_lu == NULL) {
                printf("Ошибка выделения памяти при эксперименте с несколькими правыми частями.\n");
                free(x_exact);
                free(b);
                free(x_gauss);
                free(x_lu);
                matrix_destroy(a);
                free(l);
                free(u);
                return;
            }

            gauss_task.n = n;
            gauss_task.a = a->data;
            gauss_task.b = b;
            gauss_task.x = x_gauss;
            gauss_task.mode = GAUSS_MODE_PARTIAL_PIVOT;
            gauss_task.status = GAUSS_OK;

            elapsed = measure_seconds(run_gauss_task, &gauss_task);
            gauss_total += elapsed;
            if (gauss_task.status != GAUSS_OK) {
                gauss_status = gauss_task.status;
            }

            lu_task.n = n;
            lu_task.l = l;
            lu_task.u = u;
            lu_task.b = b;
            lu_task.x = x_lu;
            lu_task.status = LU_OK;

            elapsed = measure_seconds(run_lu_factored_task, &lu_task);
            lu_total += elapsed;
            if (lu_task.status != LU_OK) {
                lu_status = lu_task.status;
            }

            free(x_exact);
            free(b);
            free(x_gauss);
            free(x_lu);
        }

        print_rhs_table_row(n, k, "gauss_partial_pivot", 0.0, gauss_total, gauss_total,
                    gauss_total / (double)k, gauss_status);

        printf("%zu,%zu,lu_factored,%d,%.6e,%.6e,%.6e,%.6e\n",
               n, k, lu_status, decompose_time, lu_total, decompose_time + lu_total,
               (decompose_time + lu_total) / (double)k);

        matrix_destroy(a);
        free(l);
        free(u);
    }
}

static void run_hilbert_experiment(void) {
    const size_t sizes[] = {5, 10, 15};
    const size_t count = sizeof(sizes) / sizeof(sizes[0]);
    size_t index;

    printf("n,method,status,time_sec,residual,relative_error\n");

    for (index = 0; index < count; ++index) {
        size_t n = sizes[index];
        Matrix *a;
        double *x_exact;
        double *b;
        double *x;

        a = generate_hilbert_matrix(n);
        x_exact = generate_constant_vector(n, 1.0);
        b = (a != NULL && x_exact != NULL) ? matrix_vector_multiply(a, x_exact) : NULL;
        x = allocate_vector(n);

        if (a == NULL || x_exact == NULL || b == NULL || x == NULL) {
            printf("%zu,error,%d,0.000000e+00,-1.000000e+00,-1.000000e+00\n", n, -1);
            matrix_destroy(a);
            free(x_exact);
            free(b);
            free(x);
            continue;
        }

        {
            GaussTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.mode = GAUSS_MODE_NO_PIVOT;
            task.status = GAUSS_OK;

            elapsed = measure_seconds(run_gauss_task, &task);
            residual = (task.status == GAUSS_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == GAUSS_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,gauss_no_pivot,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        {
            GaussTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.mode = GAUSS_MODE_PARTIAL_PIVOT;
            task.status = GAUSS_OK;

            elapsed = measure_seconds(run_gauss_task, &task);
            residual = (task.status == GAUSS_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == GAUSS_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,gauss_partial_pivot,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        {
            LUTask task;
            double elapsed;
            double residual;
            double error;

            task.n = n;
            task.a = a->data;
            task.b = b;
            task.x = x;
            task.status = LU_OK;

            elapsed = measure_seconds(run_lu_task, &task);
            residual = (task.status == LU_OK) ? residual_norm(a, x, b) : -1.0;
            error = (task.status == LU_OK) ? relative_error(n, x_exact, x) : -1.0;
            printf("%zu,lu,%d,%.6e,%.6e,%.6e\n", n, task.status, elapsed, residual, error);
        }

        matrix_destroy(a);
        free(x_exact);
        free(b);
        free(x);
    }
}

int main(void) {
    int choice;

    generate_matrix_seed(12345u);

    do {
        print_main_menu();
        if (scanf("%d", &choice) != 1) {
            printf("Ошибка ввода.\n");
            break;
        }

        if (choice == 1) {
            run_manual_case();
        } else if (choice == 2) {
            run_random_timing_experiment();
        } else if (choice == 3) {
            run_multi_rhs_experiment();
        } else if (choice == 4) {
            run_hilbert_experiment();
        } else if (choice != 0) {
            printf("Неизвестная команда.\n");
        }
    } while (choice != 0);

    return 0;
}
