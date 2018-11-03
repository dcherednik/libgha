#include <include/libgha.h>

#include "common.h"

void usage(const char* selfname) {
	fprintf(stderr, "GHA simple utility, usage: %s <FILE> <OFFSET> <LEN> [<EXPECTED_FREQ> <EXPECTED_PHASE> <EXPECTED_MAGNITUDE>]\n", selfname);
	fprintf(stderr, "FILE - input pcm file (raw pcm, mono, 24 bit, 44100 kHz\n");
	fprintf(stderr, "OFFSET - offset in input file, samples\n");
	fprintf(stderr, "LEN - len to process, samples\n");
	fprintf(stderr, "Following options are optional used by test framework:\n");
	fprintf(stderr, "EXPECTED_FREQ - expected angular frequency\n");
	fprintf(stderr, "EXPECTED_PHASE - expected phase (0 - 2Pi)\n");
	fprintf(stderr, "EXPECTED_MAGNITUDE - expected magnitude (0 - 1)\n");
}


int main(int argc, char** argv) {
	if (argc != 4 && argc != 7) {
		usage(argv[0]);
		return 1;
	}

	long long len = atoll(argv[3]);

	gha_ctx_t ctx;

	float* buf = malloc(len * sizeof(float));
	if (!buf)
		abort();

	if (load_file(argv[1], len, atoi(argv[2]), 24, buf)) {
		free(buf);
		return 1;
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
		if (fabs(freq - res.frequency) > 0.001 || compare_phase(phase, res.phase, 0.001) || fabs(magn - res.magnitude) > 0.001)
			return 1;
		return 0;
	} else {
	    fprintf(stderr, "Result: freq: %f, phase: %f, magn: %f\n", res.frequency, res.phase, res.magnitude);
	}

	return 0;
}
