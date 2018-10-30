#include <include/libgha.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void usage(const char* selfname) {
	fprintf(stderr, "GHA utility, usage: %s <FILE> <OFFSET> <LEN> [<EXPECTED_FREQ> <EXPECTED_PHASE> <EXPECTED_MAGNITUDE>]\n", selfname);
	fprintf(stderr, "FILE - input pcm file (raw pcm, mono, 24 bit, 44100 kHz\n");
	fprintf(stderr, "OFFSET - offset in input file, samples\n");
	fprintf(stderr, "LEN - len to process, samples\n");
	fprintf(stderr, "Following options are optional used by test framework:\n");
	fprintf(stderr, "EXPECTED_FREQ - expected angular frequency\n");
	fprintf(stderr, "EXPECTED_PHASE - expected phase (0 - 2Pi)\n");
	fprintf(stderr, "EXPECTED_MAGNITUDE - expected magnitude (0 - 1)\n");
}

static int compare_phase(float a, float b, float delta) {
	if (abs(a - b) < delta)
		return 0;
	a = fmod(a + M_PI, 2 * M_PI);
	b = fmod(b + M_PI, 2 * M_PI);
//	fprintf(stderr, "%f %f  %f\n", a, b, delta);
	if (abs(a - b) < delta)
		return 0;
	return 1;
}

int main(int argc, char** argv) {
	if (argc != 4 && argc != 7) {
		usage(argv[0]);
		return 1;
	}
	union {
	    char b[4];
	    int32_t i;
	} sample;

	long long i = 0;

	const char* input_file = argv[1];
	long long offset = atoll(argv[2]);
	long long len = atoll(argv[3]);
	FILE* file = fopen(input_file, "r");
	long long file_size;

	gha_ctx* ctx;

	fseek(file, 0 , SEEK_END);
	file_size = ftell(file);
	rewind(file);

	if (file_size < (offset + len) * 3) {
		return 1;
	}

	if (fseek(file, offset * 3, SEEK_SET)) {
		return 1;
	}

	float* buf = malloc(len * sizeof(float));
	if (!buf)
		abort();


	memset(sample.b, 0, sizeof(sample.b)/sizeof(sample.b[0]));
	// sample.b + 1 - sign bit at right place
	while ((fread(sample.b + 1, 1, 3, file) == 3) && i < len) {
		buf[i] = (double)(sample.i) / (double)(1u<<31);
                i++;
	}

	ctx = gha_create_ctx(len);
	if (!ctx) {
		fprintf(stderr, "Unable to create gha ctx\n");
		free(buf);
		return 1;
	}

	struct gha_info res;
	gha_analyze_one(buf, &res, ctx);

	gha_free_ctx(ctx);
	free(buf);

	if (argc == 7) {
		double freq = atof(argv[4]);
		double phase = atof(argv[5]);
		double magn = atof(argv[6]);
		if (fabs(freq - res.freq) > 0.001 || compare_phase(phase, res.phase, 0.001) || fabs(magn - res.magnitude) > 0.001)
			return 1;
		return 0;
	} else {
	    fprintf(stderr, "Result: freq: %f, phase: %f, magn: %f\n", res.freq, res.phase, res.magnitude); 
	}

	return 0;
}
