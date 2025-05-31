#include "sle.h"
#include "gha_math.h"

#include <include/libgha.h> 

#include <tools/kiss_fftr.h>

#include <stdlib.h>

#ifdef LIBGHA_PLATFORM_WINDOWS

#define alloca _alloca

#else

#ifdef LIBGHA_HAVE_ALLOCA_H

#include <alloca.h>

#endif

#endif

/*
 * Ref: http://www.apsipa.org/proceedings_2009/pdf/WA-L3-3.pdf
 */

struct gha_ctx {
	size_t size;
	size_t max_loops;
	kiss_fftr_cfg fftr;
	kiss_fftr_cfg fftr_inv;

	kiss_fft_cpx* fft_out;
	FLOAT* freq;
	FLOAT* window;

	FLOAT* tmp_buf;
	FLOAT max_magnitude;
	int upsample;
};

static void gha_init_window(gha_ctx_t ctx)
{
	size_t i;
	const size_t n = ctx->size + 1;
	const size_t half = ctx->size / 2;

	for (i = 0; i < half; i++) {
		ctx->window[i] = sinf(M_PI * (i + 1) / n);
		ctx->window[i] *= ctx->window[i];
	}

	for (i = half; i < ctx->size; i++) {
		ctx->window[i] = ctx->window[ctx->size - 1 - i];
	}
}

gha_ctx_t gha_create_ctx(size_t size)
{
	gha_ctx_t ctx = malloc(sizeof(struct gha_ctx));
	if (!ctx)
		return NULL;

	ctx->size = size;
	ctx->max_loops = 7;
	ctx->max_magnitude = 1;
	ctx->upsample = 0;

	ctx->fftr = kiss_fftr_alloc(size, 0, NULL, NULL);
	if (!ctx->fftr)
		goto exit_free_gha_ctx;

	ctx->fftr_inv = kiss_fftr_alloc(size * 2, 1, NULL, NULL);
	if (!ctx->fftr_inv)
		goto exit_free_fftr_ctx;

	ctx->freq = malloc(sizeof(FLOAT) * size);
	if (!ctx->freq)
		goto exit_free_fftr_inv_ctx;

	ctx->window = malloc(sizeof(FLOAT) * size);
	if (!ctx->window)
		goto exit_free_freq;

	ctx->tmp_buf = malloc(sizeof(FLOAT) * size);
	if (!ctx->tmp_buf)
		goto exit_free_window;

	ctx->fft_out = calloc(size + 1, sizeof(kiss_fft_cpx));
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
exit_free_fftr_inv_ctx:
	kiss_fftr_free(ctx->fftr_inv);
exit_free_fftr_ctx:
	kiss_fftr_free(ctx->fftr);
exit_free_gha_ctx:
	free(ctx);
	return NULL;
}

void gha_set_max_loops(gha_ctx_t ctx, size_t max_loops)
{
	ctx->max_loops = max_loops;
}

void gha_set_max_magnitude(gha_ctx_t ctx, FLOAT magnitude)
{
	ctx->max_magnitude = magnitude;
}

void gha_set_upsample(gha_ctx_t ctx, int enable)
{
	ctx->upsample = enable;
}

void gha_free_ctx(gha_ctx_t ctx)
{
	free(ctx->fft_out);
	free(ctx->tmp_buf);
	free(ctx->window);
	free(ctx->freq);
	kiss_fftr_free(ctx->fftr_inv);
	kiss_fftr_free(ctx->fftr);
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

static void resample_fft(gha_ctx_t ctx, FLOAT* resampled)
{
	size_t i = 0;
	kiss_fftri(ctx->fftr_inv, ctx->fft_out, resampled);
	for (i = 0; i < ctx->size * 2; i++) {
		resampled[i] /= (FLOAT)ctx->size;
	}
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

	const size_t MAX_LOOPS = 7;
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

int gha_adjust_info_newton_md(const FLOAT* pcm, struct gha_info* info, size_t dim, gha_ctx_t ctx, size_t sz)
{
	size_t loop;
	size_t i, j, k, n;

	size_t Msz = dim * 3 * (dim * 3 + 1) * sizeof(double);
	double* M = alloca(Msz);

	size_t fx0sz = dim * 3 * sizeof(double);
	double* fx0 = alloca(fx0sz);

	double* BA = alloca(sizeof(double) * dim * sz);
	double* Bw = alloca(sizeof(double) * dim * sz);
	double* Bp = alloca(sizeof(double) * dim * sz);
	double* BAw = alloca(sizeof(double) * dim * sz);
	double* BAp = alloca(sizeof(double) * dim * sz);
	double* Bww = alloca(sizeof(double) * dim * sz);
	double* Bwp = alloca(sizeof(double) * dim * sz);
	// double here breaks precision if we have only float in work buffer
	FLOAT* Bpp = alloca(sizeof(FLOAT) * dim * sz);

	for (loop = 0; loop < ctx->max_loops; loop++) {
		memcpy(ctx->tmp_buf, pcm, sz * sizeof(FLOAT));

		for (k = 0; k < dim; k++) {
			double* ba = BA + (k * sz);
			double* bw = Bw + (k * sz);
			double* bp = Bp + (k * sz);
			double* baw = BAw + (k * sz);
			double* bap = BAp + (k * sz);
			double* bww = Bww + (k * sz);
			double* bwp = Bwp + (k * sz);
			FLOAT* bpp = Bpp + (k * sz);

			for (n = 0; n < sz; n++) {
				double Ak = (info+k)->magnitude;
				float t = (info+k)->frequency * n + (info+k)->phase;
				FLOAT s = sinf(t);
				FLOAT c = cosf(t);

				ctx->tmp_buf[n] -= (info+k)->magnitude * s;

				ba[n] = -s;
				bw[n] = -Ak * n * c;
				bp[n] = -Ak * c;

				baw[n] = -n * c;
				bap[n] = -c;
				bww[n] = Ak * n * n * s;
				bwp[n] = Ak * n * s;
				bpp[n] = Ak * s;
			}
		}

		memset(M, '\0', Msz);
		for (i = 0; i < dim; i++) {
			double* m0 = M + (dim * 3 + 1) * (i + dim * 0);
			double* m1 = M + (dim * 3 + 1) * (i + dim * 1);
			double* m2 = M + (dim * 3 + 1) * (i + dim * 2);

			double* ba = BA + (i * sz);
			double* bw = Bw + (i * sz);
			double* bp = Bp + (i * sz);
			double* baw = BAw + (i * sz);
			double* bap = BAp + (i * sz);
			double* bww = Bww + (i * sz);
			double* bwp = Bwp + (i * sz);
			FLOAT* bpp = Bpp + (i * sz);

			for (j = 0; j < dim; j++) {
				if (i == j) {
					for (n = 0; n < sz; n++) {
						m0[j + dim * 0] += ba[n] * ba[n];
						m0[j + dim * 1] += ctx->tmp_buf[n] * baw[n] + ba[n] * bw[n];
						m0[j + dim * 2] += ctx->tmp_buf[n] * bap[n] + ba[n] * bp[n];

						m1[j + dim * 1] += ctx->tmp_buf[n] * bww[n] + bw[n] * bw[n];
						m1[j + dim * 2] += ctx->tmp_buf[n] * bwp[n] + bw[n] * bp[n];

						m2[j + dim * 2] += ctx->tmp_buf[n] * bpp[n] + bp[n] * bp[n];
					}
				} else {
					for (n = 0; n < sz; n++) {
						double* baj = BA + (j * sz);
						double* bpj = Bp + (j * sz);
						double* bwj = Bw + (j * sz);
						m0[j + dim * 0] += ba[n] * baj[n];
						m0[j + dim * 1] += ba[n] * bwj[n];
						m0[j + dim * 2] += ba[n] * bpj[n];

						m1[j + dim * 1] += bw[n] * bwj[n];
						m1[j + dim * 2] += bw[n] * bpj[n];

						m2[j + dim * 2] += bp[n] * bpj[n];
					}
				}
				m0[j + dim * 0] *= 2;
				m0[j + dim * 1] *= 2;
				m0[j + dim * 2] *= 2;

				m1[j + dim * 1] *= 2;
				m1[j + dim * 2] *= 2;

				m2[j + dim * 2] *= 2;


				m0[j + dim * 1] = m1[j + dim * 0];
				m0[j + dim * 2] = m2[j + dim * 0];
				m1[j + dim * 2] = m2[j + dim * 1];
			}
		}

		for (k = 0; k < dim; k++) {
			double* m0 = M + (dim * 3 + 1) * (k + dim * 0);
			double* m1 = M + (dim * 3 + 1) * (k + dim * 1);
			double* m2 = M + (dim * 3 + 1) * (k + dim * 2);
			double* ba = BA + (k * sz); 
			double* bw = Bw + (k * sz); 
			double* bp = Bp + (k * sz); 
			for (n = 0; n < sz; n++) {
				m0[dim * 3] += ctx->tmp_buf[n] * (FLOAT)ba[n];
				m1[dim * 3] += ctx->tmp_buf[n] * (FLOAT)bw[n];
				m2[dim * 3] += ctx->tmp_buf[n] * (FLOAT)bp[n];
			}
			m0[dim * 3] *= 2;
			m1[dim * 3] *= 2;
			m2[dim * 3] *= 2;
		}

		memset(fx0, '\0', fx0sz);
		if(sle_solve(M, dim * 3, fx0)) {
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

			if ((info+k)->magnitude > ctx->max_magnitude) {
				//TODO: ???
				(info+k)->magnitude = ctx->max_magnitude * 0.5;
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

	FLOAT* resampled;

	if (ctx->upsample > 0) {
		resampled = alloca(sizeof(FLOAT) * ctx->size * 2);

		resample_fft(ctx, resampled);

		gha_search_omega_newton(resampled, bin, ctx->size * 2, info);
		info->frequency *= 2.0;
	} else {
		gha_search_omega_newton(ctx->tmp_buf, bin, ctx->size, info);
	}

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
}

void gha_extract_many_simple(FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx)
{
	int i;
	for (i = 0; i < k; i++) {
		gha_extract_one(pcm, info + i, ctx);
	}
}

int gha_adjust_info(const FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx, resuidal_cb_t cb, void* user_ctx, size_t size_limit)
{
	size_t actual_size = ctx->size;
	int rv;

	if (size_limit && size_limit < ctx->size) {
		actual_size = size_limit;
	}

	rv = gha_adjust_info_newton_md(pcm, info, k, ctx, actual_size);
	if (cb && rv != -1)
		cb(ctx->tmp_buf, ctx->size, user_ctx);

	return rv;
}

const FLOAT* gha_get_analyzed(gha_ctx_t ctx)
{
	return ctx->tmp_buf;
}
