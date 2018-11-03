#include <include/libgha.h>

#include "common.h"

void usage(const char* selfname) {
	fprintf(stderr, "GHA dtmf example utility, usage: %s <FILE> <OFFSET> <LEN> [<EXPECTED_FREQ_LOW> <EXPECTED_MAGNITUDE_LOW>, <EXPECTED_FREQ_HIGH>, <EXPECTED_MAGNITUDE_HIGH>]\n", selfname);
	fprintf(stderr, "FILE - input pcm file (raw pcm, mono, 8 bit, 8000 kHz\n");
	fprintf(stderr, "OFFSET - offset in input file, samples\n");
	fprintf(stderr, "LEN - len to process, samples\n");
	fprintf(stderr, "Following options are optional used by test framework:\n");
	fprintf(stderr, "EXPECTED_FREQ_(LOW|HIGH) - expected angular frequency\n");
	fprintf(stderr, "EXPECTED_MAGNITUDE_(LOW|HIGH) - expected magnitude (0 - 1)\n");
}

int main(int argc, char** argv) {
	if (argc != 4 && argc != 8) {
		usage(argv[0]);
		return 1;
	}

	long long len = atoll(argv[3]);

	gha_ctx_t ctx;

	float* buf = malloc(len * sizeof(float));
	if (!buf)
		abort();

	if (load_file(argv[1], len, atoi(argv[2]), 8, buf)) {
		fprintf(stderr, "Unable to load input file\n");
		free(buf);
		return 1;
	}


	ctx = gha_create_ctx(len);
	if (!ctx) {
		fprintf(stderr, "Unable to create gha ctx\n");
		free(buf);
		return 1;
	}

	struct gha_info res[2];
	gha_extract_one(buf, &res[0], ctx);
	gha_extract_one(buf, &res[1], ctx);

	if (res[0].frequency > res[1].frequency) {
		struct gha_info tmp;
		memcpy(&tmp, &res[0], sizeof(struct gha_info));
		memcpy(&res[0], &res[1], sizeof(struct gha_info));
		memcpy(&res[1], &tmp, sizeof(struct gha_info));
	}

	gha_free_ctx(ctx);
	free(buf);

	if (argc == 8) {
		double freq[2];
		double magn[2];
		freq[0] = atof(argv[4]);
		magn[0] = atof(argv[5]);
		freq[1] = atof(argv[6]);
		magn[1] = atof(argv[7]);
		if (fabs(freq[0] - res[0].frequency) > 0.001 || fabs(magn[0] - res[0].magnitude) > 0.001 ||
			fabs(freq[1] - res[1].frequency) > 0.001 || fabs(magn[1] - res[1].magnitude) > 0.001)
			return 1;
		return 0;


	} else {
	    fprintf(stderr, "dtmf result: low freq: %f, magn: %f, high freq: %f, magn: %f\n",
			res[0].frequency, res[0].magnitude, res[1].frequency, res[1].magnitude);
	}
}
