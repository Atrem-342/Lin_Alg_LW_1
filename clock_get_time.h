#ifndef CLOCK_GET_TIME_H
#define CLOCK_GET_TIME_H

#include <stddef.h>

typedef void (*TimeFunction)(void *arg);

double now_seconds(void);
double measure_seconds(TimeFunction func, void *arg);

#endif
