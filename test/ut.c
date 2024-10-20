#include <3rd/fctx/fct.h>

#include <sle.h>
#include "common.h"

static double eq_matrix_1[3][4] =
	{{ 2,	 1,	-1,	 8},
	{-3,	-1,	 2,	-11},
	{-2,	 1,	 2,	-3}};
static double eq_result_1[3] =
	{2.0, 3.0, -1.0};

static double eq_matrix_2[2][3] =
	{{ 2,	 1,	1},
	{4,	2,	2}};

#define SAMPLERATE 44100

static void gen(float f, float a, float* out, size_t sz)
{
	float freq = f / SAMPLERATE;
	size_t i;
	for (i = 0; i < sz; i++) {
		out[i] += sin(freq * (float)(i * 2 * M_PI)) * a;
	}
}


FCT_BGN()
{
	FCT_SUITE_BGN(internal)
	{
		FCT_TEST_BGN(sle_ok_1)
		{
			double* result = calloc(3, sizeof(double));;
			int i, rv;
			rv = sle_solve(eq_matrix_1[0], 3, result);
			fct_chk_eq_int(rv, 0);
			for (i = 0; i < 3; i++)
				fct_chk_eq_dbl(result[i], eq_result_1[i]);
			free(result);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(sle_not_ok_2)
		{
			double* result = calloc(2, sizeof(double));;
			int rv;
			rv = sle_solve(eq_matrix_2[0], 2, result);
			fct_chk_eq_int(rv, -1);
			free(result);
		}
		FCT_TEST_END();
	}
	FCT_SUITE_END();

	FCT_SUITE_BGN(gha)
		FCT_TEST_BGN(one_tone_11025_a1)
		{
			float  buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res;
			ctx = gha_create_ctx(128);

			gen(11025.0, 1, buf, 128);

			gha_analyze_one(buf, &res, ctx);

			fprintf(stderr, "Result: freq: %.10f, phase: %f, magn: %f\n", res.frequency, res.phase, res.magnitude);
			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 1.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 0);
			fprintf(stderr, "Result: freq: %.10f, phase: %f, magn: %f\n", res.frequency, res.phase, res.magnitude);
			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 1.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			res.magnitude = 0.95; /* add some error */
			/* asume this routine will be able to fix error using just part of buffer */
			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 32);
			fprintf(stderr, "Result: freq: %.10f, phase: %f, magn: %f\n", res.frequency, res.phase, res.magnitude);
			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 1.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);


			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(one_tone_11025_a32768)
		{
			float  buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res;
			ctx = gha_create_ctx(128);
			gha_set_max_magnitude(ctx, 32768);

			gen(11025.0, 32768, buf, 128);

			gha_analyze_one(buf, &res, ctx);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 32768.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 0);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 32768.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			res.magnitude = 32760.0;

			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 32);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 32768.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(one_tone_11025_a32768_adjust)
		{
			float  buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res;
			res.magnitude = 1000;
			res.phase = 0.0;
			res.frequency = 1.5;

			ctx = gha_create_ctx(128);
			gha_set_max_loops(ctx, 128);
			gha_set_max_magnitude(ctx, 32768);

			gen(11025.0, 32768, buf, 128);

			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 0);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 32768.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(two_tones_5000_11025_a32768_adjust)
		{
			float  buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res[2];
			res[0].magnitude = 16000;
			res[0].phase = 0.0;
			res[0].frequency = 1.5;

			res[1].magnitude = 16000;
			res[1].phase = 0.0;
			res[1].frequency = 0.75;

			ctx = gha_create_ctx(128);
			gha_set_max_loops(ctx, 128);
			gha_set_max_magnitude(ctx, 32768);

			gen(11025.0, 32768, buf, 128);
			gen(5000.0, 32768, buf, 128);

			gha_adjust_info(buf, res, 2, ctx, NULL, NULL, 0);

			fct_chk_eq_int(0, compare_phase(0, res[0].phase, 0.01));
			UT_CHECK_EQ_FLOAT(round(res[0].magnitude), 32768.0);
			UT_CHECK_EQ_FLOAT(res[0].frequency, 1.5707963705);

			fct_chk_eq_int(0, compare_phase(0, res[1].phase, 0.01));
			UT_CHECK_EQ_FLOAT(round(res[1].magnitude), 32768.0);
			UT_CHECK_EQ_FLOAT(res[1].frequency, 0.712379);

			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(two_tones_5512hz5_11025_a32768_adjust)
		{
			float  buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res[2];
			res[0].magnitude = 16000;
			res[0].phase = 0.0;
			res[0].frequency = 1.5707;

			res[1].magnitude = 16000;
			res[1].phase = 0.0;
			res[1].frequency = 0.75;

			ctx = gha_create_ctx(128);
			gha_set_max_loops(ctx, 128);
			gha_set_max_magnitude(ctx, 32768);

			gen(11025.0, 16384, buf, 128);
			gen(5512.5, 32768, buf, 128);

			gha_adjust_info(buf, res, 2, ctx, NULL, NULL, 0);

			fct_chk_eq_int(0, compare_phase(0, res[0].phase, 0.01));
			UT_CHECK_EQ_FLOAT(round(res[0].magnitude), 16384.0);
			UT_CHECK_EQ_FLOAT(res[0].frequency, 1.5707963705);

			fct_chk_eq_int(0, compare_phase(0, res[1].phase, 0.01));
			UT_CHECK_EQ_FLOAT(round(res[1].magnitude), 32768.0);
			UT_CHECK_EQ_FLOAT(res[1].frequency, 0.785398);

			fprintf(stderr, "Result: freq: %.10f, phase: %f, magn: %f\n", res[0].frequency, res[0].phase, res[0].magnitude);
			fprintf(stderr, "Result: freq: %.10f, phase: %f, magn: %f\n", res[1].frequency, res[1].phase, res[1].magnitude);

			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(one_tone_11025_a32768_partial_frame)
		{
			float buf[128] = {0};
			gha_ctx_t ctx;
			struct gha_info res;
			ctx = gha_create_ctx(128);
			gha_set_max_magnitude(ctx, 32768);
			gha_set_max_loops(ctx, 14);

			gen(11025.0, 32768, buf, 128);

			memset(&buf[96], '\0', 32*sizeof(float));

			gha_analyze_one(buf, &res, ctx);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 24576.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_adjust_info(buf, &res, 1, ctx, NULL, NULL, 64);

			fct_chk_eq_int(0, compare_phase(0, res.phase, 0.01));
			UT_CHECK_EQ_FLOAT(res.magnitude, 32768.0);
			UT_CHECK_EQ_FLOAT(res.frequency, 1.5707963705);

			gha_free_ctx(ctx);
		}
		FCT_TEST_END();

	FCT_SUITE_END();
}
FCT_END();
