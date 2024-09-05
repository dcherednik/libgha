#ifndef LIBGHA_H
#define LIBGHA_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#ifdef GHA_USE_DOUBLE_API
#	define kiss_fft_scalar double
#	define FLOAT double
#else
#	define FLOAT float
#endif

#include <stddef.h>

typedef struct gha_ctx *gha_ctx_t;
typedef void (*resuidal_cb_t)(FLOAT* resuidal, size_t size, void* user_ctx);

struct gha_info {
	FLOAT frequency;
	FLOAT phase;
	FLOAT magnitude;
};

/*
 * Create context to perform GHA, size is number of samples provided to analyze.
 * Size must be even
 *
 * Returns null in case of fail.
 *
 */
gha_ctx_t gha_create_ctx(size_t size);

void gha_set_max_loops(gha_ctx_t ctx, size_t max_loops);
void gha_set_max_magnitude(gha_ctx_t ctx, FLOAT magnitude);

/*
 * Free GHA context
 */
void gha_free_ctx(gha_ctx_t ctx);

/*
 * Performs one GHA step for given PCM signal,
 * the result will be writen in to given gha_info structure
 *
 * Complexity: O(n * log(n)),
 * where n is number of samples to anayze
 *
 */
void gha_analyze_one(const FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx);

/*
 * Performs one GHA step and extracts analysed harmonic from given PCM signal
 * the result will be writen in to given gha_info structure
 *
 * Complexity: O(n * log(n)),
 * where n is number of samples to anayze
 *
 */
void gha_extract_one(FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx);

/*
 * Performs k GHA steps and extracts analysed harmonics from given PCM signal.
 *
 * The result will be writen in to given gha_info structures
 * corresponding ammount of memory should be allocated (k * sizeof(struct gha_info))
 *
 * Effectively this function is equivalent of calling gha_extract_one k times
 *
 * Complexity: O(n * log(n) * k),
 * where n is number of samples to anayze, k is number of harmonics to extract
 *
 */
void gha_extract_many_simple(FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx);

/*
 * Performs multidimensional optimization of extracted harmonics.
 *
 * Given gha_info will be adjusted to minimize resuidal level.
 * Newton multidimensional optimization method is used.
 *
 * Complexity: O(k^2 * n + k^3)
 * where n is number of samples to anayze, k is number of harmonics to extract
 *
 */
int gha_adjust_info(const FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx, resuidal_cb_t cb, void* user_ctx);

const FLOAT* gha_get_analyzed(gha_ctx_t ctx);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
