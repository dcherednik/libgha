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

static void calc_resuidal(FLOAT* resuidal, size_t size, void* user_ctx)
{
	int i;
	double s = 0.0;
	FLOAT* result = (FLOAT*)user_ctx;
	for (i = 0; i < size; i++) {
		s += resuidal[i] * resuidal[i];
	}
	*result = sqrt(s / size);
}

int main(int argc, char** argv) {
	if (argc != 4 && argc != 8) {
		usage(argv[0]);
		return 1;
	}

	long long len = atoll(argv[3]);

	gha_ctx_t ctx;

	FLOAT* buf = malloc(len * sizeof(FLOAT));
	FLOAT* buf2 = malloc(len * sizeof(FLOAT));
	if (!buf || !buf2)
		abort();

	if (load_file(argv[1], len, atoi(argv[2]), 8, buf)) {
		fprintf(stderr, "Unable to load input file\n");
		free(buf);
		return 1;
	}

	//Make copy of data to adjust extracted params
	memcpy(buf2, buf, sizeof(FLOAT) * len);
	ctx = gha_create_ctx(len);

	FLOAT resuidal;
	gha_set_user_resuidal_cb(&calc_resuidal, &resuidal, ctx);
	if (!ctx) {
		fprintf(stderr, "Unable to create gha ctx\n");
		free(buf);
		return 1;
	}

	struct gha_info res[2];
	gha_extract_many_simple(buf, &res[0], 2, ctx);

	FLOAT resuidal_1 = resuidal;

	if (res[0].frequency > res[1].frequency) {
		struct gha_info tmp;
		memcpy(&tmp, &res[0], sizeof(struct gha_info));
		memcpy(&res[0], &res[1], sizeof(struct gha_info));
		memcpy(&res[1], &tmp, sizeof(struct gha_info));
	}

	gha_adjust_info(buf2, res, 2, ctx);

	if (resuidal > resuidal_1) {
		fprintf(stderr, "gha_adjust_info wrong result\n");
		return 1;
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
	    return 0;
	}
}
