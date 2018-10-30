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
};

static void gha_init_window(gha_ctx* ctx)
{
	size_t i;
	size_t n = ctx->size + 1;
	for (i = 0; i < ctx->size; i++) {
		ctx->window[i] = sin(M_PI * (i + 1) / n);
	}
}

gha_ctx* gha_create_ctx(size_t size)
{
	gha_ctx* ctx = malloc(sizeof(struct gha_ctx));
	if (!ctx)
		return NULL;

	ctx->size = size;

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

void gha_free_ctx(gha_ctx* ctx)
{
	free(ctx->fft_out);
	free(ctx->tmp_buf);
	free(ctx->window);
	free(ctx->freq);
	kiss_fft_free(ctx->fftr);
	free(ctx);
}

static size_t gha_estimate_bin(gha_ctx* ctx)
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

		for (n = 0; n < size; n++) {
			double c = pcm[n] * cos(omega_rad * n);
			double s = pcm[n] * sin(omega_rad * n);
			double tc, ts;
			Xr += c;
			Xi += s;
			tc = n * c;
			ts = n * s;
			dXr -= ts;
			dXi += tc;
			ddXr -= n * tc;
			ddXs -= n * ts;
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
		    result->freq = omega_rad;
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

void gha_analyze_one(const FLOAT* pcm, struct gha_info* info, gha_ctx* ctx)
{
	int i = 0;
	int bin = 0;

	for (i = 0; i < ctx->size; i++)
		ctx->tmp_buf[i] = pcm[i] * ctx->window[i];

	kiss_fftr(ctx->fftr, ctx->tmp_buf, ctx->fft_out);

	bin = gha_estimate_bin(ctx);

	gha_search_omega_newton(ctx->tmp_buf, bin, ctx->size, info);
	gha_generate_sine(ctx->tmp_buf, ctx->size, info->freq, info->phase);
	gha_estimate_magnitude(pcm, ctx->tmp_buf, ctx->size, info);
}
