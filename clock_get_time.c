#define _POSIX_C_SOURCE 200809L

#include "clock_get_time.h"

#include <time.h>

double now_seconds(void) {
    struct timespec ts;

    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0;
}

double measure_seconds(TimeFunction func, void *arg) {
    double start;
    double end;

    if (func == NULL) {
        return 0.0;
    }

    start = now_seconds();
    func(arg);
    end = now_seconds();

    return end - start;
}
