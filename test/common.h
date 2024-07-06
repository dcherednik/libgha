#include <include/libgha.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static inline int compare_phase(FLOAT a, FLOAT b, FLOAT delta) {
	if (fabs(a - b) < delta)
		return 0;
	a = fmod(a + M_PI, 2 * M_PI);
	b = fmod(b + M_PI, 2 * M_PI);
//	fprintf(stderr, "%f %f  %f\n", a, b, delta);
	if (fabs(a - b) < delta)
		return 0;
	return -1;
}

static inline int load_file(const char* name, size_t len, size_t offset, size_t bits, FLOAT* buf)
{
	union {
	    char b[4];
	    int32_t i;
	} sample;
	int bytes_per_sample = 0;

	size_t file_size, i;

	FILE* file = fopen(name, "r");

	if (!file)
		return -1;

	switch (bits) {
		case 8:
			bytes_per_sample = 1;
			break;
		case 24:
			bytes_per_sample = 3;
			break;
		default:
			fprintf(stderr, "Unsupported sample format\n");
			return -1;
	}

	fseek(file, 0, SEEK_END);
	file_size = ftell(file);
	rewind(file);

	if (file_size < (offset + len) * bytes_per_sample)
		return -1;

	if (fseek(file, offset * bytes_per_sample, SEEK_SET)) {
		return -1;
	}

	memset(sample.b, 0, sizeof(sample.b)/sizeof(sample.b[0]));

	i = 0;
	char* p = sample.b + 4 - bytes_per_sample;
	while ((fread(p, 1, bytes_per_sample, file) == bytes_per_sample) && i < len) {
		buf[i] = (double)(sample.i) / (double)(1u<<31);
                i++;
	}
	return 0;
}
