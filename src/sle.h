#ifndef SLE_H
#define SLE_H

#include <stddef.h>
/*
 * m - matrix[n][n+1]
 * n - row count
 * x - result
 * returns 0 in case of success, -1 error
 */
int sle_solve(double *a, size_t n, double *x);

#endif
