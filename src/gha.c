#include "sle.h"

#include <include/libgha.h> 

#include <kissfft/tools/kiss_fftr.h>

/*
 * Ref: http://www.apsipa.org/proceedings_2009/pdf/WA-L3-3.pdf
 */

struct gha_ctx {
	size_t size;
	kiss_fftr_cfg fftr;

	kiss_fft_cpx* fft_out;
	FLOAT* freq;
	FLOAT* window;

	FLOAT* tmp_buf;

	void (*resuidal_cb)(FLOAT* resuidal, size_t size, void* user_ctx);
	void* user_ctx;
};

static void gha_init_window(gha_ctx_t ctx)
{
	size_t i;
	size_t n = ctx->size + 1;
	for (i = 0; i < ctx->size; i++) {
		ctx->window[i] = sin(M_PI * (i + 1) / n);
	}
}

gha_ctx_t gha_create_ctx(size_t size)
{
	gha_ctx_t ctx = malloc(sizeof(struct gha_ctx));
	if (!ctx)
		return NULL;

	ctx->size = size;
	ctx->resuidal_cb = NULL;
	ctx->user_ctx = NULL;

	ctx->fftr = kiss_fftr_alloc(size, 0, NULL, NULL);
	if (!ctx->fftr)
		goto exit_free_gha_ctx;

	ctx->freq = malloc(sizeof(FLOAT) * size);
	if (!ctx->freq)
		goto exit_free_fftr_ctx;

	ctx->window = malloc(sizeof(FLOAT) * size);
	if (!ctx->window)
		goto exit_free_freq;

	ctx->tmp_buf = malloc(sizeof(FLOAT) * size);
	if (!ctx->tmp_buf)
		goto exit_free_window;

	ctx->fft_out = malloc(sizeof(kiss_fft_cpx) * (size/2 + 1));
	if (!ctx->fft_out)
		goto exit_free_tmp_buf;

	gha_init_window(ctx);

	return ctx;
exit_free_tmp_buf:
	free(ctx->tmp_buf);
exit_free_window:
	free(ctx->window);
exit_free_freq:
	free(ctx->freq);
exit_free_fftr_ctx:
	kiss_fftr_free(ctx->fftr);
exit_free_gha_ctx:
	free(ctx);
	return NULL;
}

void gha_set_user_resuidal_cb(void (*cb)(FLOAT* resuidal, size_t size, void* user_ctx), void* user_ctx, gha_ctx_t ctx)
{
	ctx->user_ctx = user_ctx;
	ctx->resuidal_cb = cb;
}

void gha_free_ctx(gha_ctx_t ctx)
{
	free(ctx->fft_out);
	free(ctx->tmp_buf);
	free(ctx->window);
	free(ctx->freq);
	kiss_fft_free(ctx->fftr);
	free(ctx);
}

static size_t gha_estimate_bin(gha_ctx_t ctx)
{
	size_t i, end;
	size_t j = 0;
	FLOAT max = 0.0;
	FLOAT tmp = 0.0;
	end = ctx->size/2 + 1;
	for (i = 0; i < end; i++) {
		tmp = ctx->fft_out[i].r * ctx->fft_out[i].r + ctx->fft_out[i].i * ctx->fft_out[i].i;
		if (tmp > max) {
			max = tmp;
			j = i;
		}
	}
	return j;
}

/*
 * Perform search of frequency using Newton's method
 * Also we calculate real and imaginary part of Fourier transform at target frequency
 * so we also calculate phase here at last iteration
 */
static void gha_search_omega_newton(const FLOAT* pcm, size_t bin, size_t size, struct gha_info* result)
{
	size_t loop;
	int n;
	double omega_rad = bin * 2 * M_PI / size;

	const size_t MAX_LOOPS = 8;
	for (loop = 0; loop <= MAX_LOOPS; loop++) {
		double Xr = 0;
		double Xi = 0;
		double dXr = 0;
		double dXi = 0;
		double ddXr = 0;
		double ddXs = 0;

		const double a = cos(omega_rad);
		const double b = sin(omega_rad);
		double c = 1.0;
		double s = 0.0;

		for (n = 0; n < size; n++) {
			double cm = pcm[n] * c;
			double sm = pcm[n] * s;
			double tc, ts;
			Xr += cm;
			Xi += sm;
			tc = n * cm;
			ts = n * sm;
			dXr -= ts;
			dXi += tc;
			ddXr -= n * tc;
			ddXs -= n * ts;

			const double new_c = a * c - b * s;
			const double new_s = b * c + a * s;
			c = new_c;
			s = new_s;


		}

		double F = Xr * dXr + Xi * dXi;
		double G2 = Xr * Xr + Xi * Xi;
		//fprintf(stderr, " %f %f \n", Xr, Xi);
		//double dXg = F;
		double dF = Xr * ddXr + dXr * dXr + Xi * ddXs + dXi * dXi;

		//double dg = F / G;
		//double ddXg = (dF * G - F * dg) / G;
		//double dw = dXg / ddXg;
		double dw = F / (dF - (F * F) / G2);
		//fprintf(stderr, "dw: %f\n", dw);

		omega_rad -= dw;

		if (omega_rad < 0)
			omega_rad *= -1;

		while (omega_rad > M_PI * 2.0)
			omega_rad -= M_PI * 2.0;

		if (omega_rad > M_PI)
			omega_rad = M_PI * 2.0 - omega_rad;

		// Last iteration
		if (loop == MAX_LOOPS) {
		    result->frequency = omega_rad;
		    //assume zero-phase sine
		    result->phase = M_PI / 2 - atan(Xi / Xr);
		    if (Xr < 0)
			    result->phase += M_PI;
		}
	}
}

static void gha_generate_sine(FLOAT* buf, size_t size, FLOAT omega, FLOAT phase)
{
	int i;
	for (i = 0; i < size; i++) {
		buf[i] = sin(omega * i + phase);
	}
}

static void gha_estimate_magnitude(const FLOAT* pcm, const FLOAT* regen, size_t size, struct gha_info* result)
{
	int i;
	double t1 = 0;
	double t2 = 0;
	for (i = 0; i < size; i++) {
		t1 += pcm[i] * regen[i];
		t2 += regen[i] * regen[i];
	}

	result->magnitude = t1 / t2;
}

int gha_adjust_info_newton_md(const FLOAT* pcm, struct gha_info* info, size_t dim, gha_ctx_t ctx)
{
	size_t loop;
	size_t i, j, k, n;

	const size_t MAX_LOOPS = 7;

	for (loop = 0; loop < MAX_LOOPS; loop++) {
		memcpy(ctx->tmp_buf, pcm, ctx->size * sizeof(FLOAT));

		// Use VLA for a while
		FLOAT BA[dim][ctx->size];
		FLOAT Bw[dim][ctx->size];
		FLOAT Bp[dim][ctx->size];
		FLOAT BAw[dim][ctx->size];
		FLOAT BAp[dim][ctx->size];
		FLOAT Bww[dim][ctx->size];
		FLOAT Bwp[dim][ctx->size];
		FLOAT Bpp[dim][ctx->size];

		for (n = 0; n < ctx->size; n++) {
			for (k = 0; k < dim; k++) {
				FLOAT Ak = (info+k)->magnitude;
				FLOAT wk = (info+k)->frequency;
				FLOAT pk = (info+k)->phase;
				FLOAT s = sin(wk * n + pk);
				FLOAT c = cos(wk * n + pk);
				ctx->tmp_buf[n] -= Ak * s;

				BA[k][n] = -s;
				Bw[k][n] = -Ak * n * c;
				Bp[k][n] = -Ak * c;

				BAw[k][n] = -n * c;
				BAp[k][n] = -c;
				Bww[k][n] = Ak * n * n * s;
				Bwp[k][n] = Ak * n * s;
				Bpp[k][n] = Ak * s;

			}
		}

		double M[dim * 3][dim * 3 + 1];
		memset(M, '\0', dim * 3 * (dim * 3 + 1) * sizeof(double));
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				for (n = 0; n < ctx->size; n++) {
					if (i == j) {
						M[i + dim * 0][j + dim * 0] += BA[i][n] * BA[i][n];
						M[i + dim * 0][j + dim * 1] += ctx->tmp_buf[n] * BAw[i][n] + BA[i][n] * Bw[i][n];
						M[i + dim * 0][j + dim * 2] += ctx->tmp_buf[n] * BAp[i][n] + BA[i][n] * Bp[i][n];

						M[i + dim * 1][j + dim * 1] += ctx->tmp_buf[n] * Bww[i][n] + Bw[i][n] * Bw[i][n];
						M[i + dim * 1][j + dim * 2] += ctx->tmp_buf[n] * Bwp[i][n] + Bw[i][n] * Bp[i][n];

						M[i + dim * 2][j + dim * 2] += ctx->tmp_buf[n] * Bpp[i][n] + Bp[i][n] * Bp[i][n];
					} else {
						M[i + dim * 0][j + dim * 0] += BA[i][n] * BA[j][n];
						M[i + dim * 0][j + dim * 1] += BA[i][n] * Bw[j][n];
						M[i + dim * 0][j + dim * 2] += BA[i][n] * Bp[j][n];

						M[i + dim * 1][j + dim * 1] += Bw[i][n] * Bw[j][n];
						M[i + dim * 1][j + dim * 2] += Bw[i][n] * Bp[j][n];

						M[i + dim * 2][j + dim * 2] += Bp[i][n] * Bp[j][n];
					}
				}
				M[i + dim * 0][j + dim * 0] *= 2;
				M[i + dim * 0][j + dim * 1] *= 2;
				M[i + dim * 0][j + dim * 2] *= 2;

				M[i + dim * 1][j + dim * 1] *= 2;
				M[i + dim * 1][j + dim * 2] *= 2;

				M[i + dim * 2][j + dim * 2] *= 2;


				M[i + dim * 0][j + dim * 1] = M[i + dim * 1][j + dim * 0];
				M[i + dim * 0][j + dim * 2] = M[i + dim * 2][j + dim * 0];
				M[i + dim * 1][j + dim * 2] = M[i + dim * 2][j + dim * 1];
			}
		}

		for (k = 0; k < dim; k++) {
			for (n = 0; n < ctx->size; n++) {
				M[k + dim * 0][dim * 3] += ctx->tmp_buf[n] * BA[k][n];
				M[k + dim * 1][dim * 3] += ctx->tmp_buf[n] * Bw[k][n];
				M[k + dim * 2][dim * 3] += ctx->tmp_buf[n] * Bp[k][n];
			}
			M[k + dim * 0][dim * 3] *= 2;
			M[k + dim * 1][dim * 3] *= 2;
			M[k + dim * 2][dim * 3] *= 2;
		}

		double fx0[dim * 3];
		memset(fx0, '\0', dim * 3 * sizeof(double));
		if(sle_solve(&M[0][0], dim * 3, fx0)) {
			return -1;
		}

		for (k = 0; k < dim; k++) {
			//fprintf(stderr, "delta1: %f\n", fx0[k + dim * 0]);
			//fprintf(stderr, "delta2: %f\n", fx0[k + dim * 1]);
			//fprintf(stderr, "delta3: %f\n", fx0[k + dim * 2]);
			(info+k)->magnitude -= (fx0[k + dim * 0] * 0.8);
			(info+k)->frequency -= (fx0[k + dim * 1] * 0.8);
			(info+k)->phase -=     (fx0[k + dim * 2] * 0.8);
		}

		for (k = 0; k < dim; k++) {
			if ((info+k)->magnitude < 0) {
				(info+k)->magnitude *= -1;
				(info+k)->phase += M_PI;
			}

			if ((info+k)->magnitude > 1) {
				//TODO: ???
				(info+k)->magnitude = 0.5;
			}
		}

		for (k = 0; k < dim; k++) {
			if ((info + k)->frequency < 0) {
				//fprintf(stderr, "negative freq\n");
				(info + k)->frequency *= -1;
				(info + k)->phase = 2 * M_PI - (info + k)->phase;
			}
			while ((info + k)->frequency > M_PI * 2.0) {
				//fprintf(stderr, "freq over\n");
				(info + k)->frequency -= M_PI * 2.0;
			}
			if ((info + k)->frequency > M_PI) {
				//fprintf(stderr, "freq ??\n");
				(info + k)->frequency = 2 * M_PI - (info + k)->frequency;
			}
		}

		for (k = 0; k < dim; k++) {
			while ((info+k)->phase > M_PI * 2.0) {
				(info+k)->phase -= M_PI * 2;
			}
			while ((info + k)->phase < 0) {
				(info+k)->phase += M_PI * 2;
			}
		}
	}
	return 0;
}

void gha_analyze_one(const FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx)
{
	int i = 0;
	int bin = 0;

	for (i = 0; i < ctx->size; i++)
		ctx->tmp_buf[i] = pcm[i] * ctx->window[i];

	kiss_fftr(ctx->fftr, ctx->tmp_buf, ctx->fft_out);

	bin = gha_estimate_bin(ctx);

	gha_search_omega_newton(ctx->tmp_buf, bin, ctx->size, info);
	gha_generate_sine(ctx->tmp_buf, ctx->size, info->frequency, info->phase);
	gha_estimate_magnitude(pcm, ctx->tmp_buf, ctx->size, info);
}

void gha_extract_one(FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx)
{
	int i;
	FLOAT magnitude;
	gha_analyze_one(pcm, info, ctx);
	magnitude = info->magnitude;

	for (i = 0; i < ctx->size; i++)
		pcm[i] -= ctx->tmp_buf[i] * magnitude;

	if (ctx->resuidal_cb)
		ctx->resuidal_cb(pcm, ctx->size, ctx->user_ctx);
}

void gha_extract_many_simple(FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx)
{
	int i;
	for (i = 0; i < k; i++) {
		gha_extract_one(pcm, info + i, ctx);
	}
}

int gha_adjust_info(const FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx)
{
	int rv = gha_adjust_info_newton_md(pcm, info, k, ctx);
	if (ctx->resuidal_cb)
		ctx->resuidal_cb(ctx->tmp_buf, ctx->size, ctx->user_ctx);

	return rv;
}
