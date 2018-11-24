#include "sle.h" 
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

/*
 * Gaussian elimination
 */
int sle_solve(double *a, size_t n, double *x) {
	int i, j, k, index;
	double max, t;
	const int col = n + 1;
	const float eps = 0.00001;

	if (n == 0) {
		return - 1;
	}
	for (k = 0; k < n; k++) {
		max = fabs(a[col*k + k]);
		index = k;
		for (i = k + 1; i < n; i++) {
			t = abs(a[col * i + k]);
			if (t > max) {
				max = t;
				index = i;
			}
		}
		if (max < eps) {
			return -1;
		}
		if (index != k) {
			for (i = 0; i < col; i++) {
				t = a[col * k + i];
				a[col * k + i] = a[col * index + i];
				a[col * index + i] = t;
			}
		}
		for (i = k; i < n; i++) {
			t = a[col * i + k];
			if (fabs(t) < eps)
				continue;
			for (j = 0; j < col; j++) {
				a[i * col + j] /= t;
			}
			if (i != k) {
				for (j = 0; j < col; j++) {
					a[i * col + j] -= a[k * col + j];
				}
			}
		}
	}

	for (k = n - 1; k >= 0; k--) {
		x[k] = a[col * k + n];
		for (i = 0; i < k; i++) {
			a[col * i + n] -= a[col * i + k] * x[k];
		}
	}
	return 0;
}
