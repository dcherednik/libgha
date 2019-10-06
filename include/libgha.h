#ifndef LIBGHA_H
#define LIBGHA_H

#define FLOAT float

#include <stddef.h>

typedef struct gha_ctx *gha_ctx_t;

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
int gha_adjust_info(const FLOAT* pcm, struct gha_info* info, size_t k, gha_ctx_t ctx);

/*
 * Set callback to perform action on resuidal pcm signal.
 *
 * user_ctx is used to pass user context in to callback
 *
 */
void gha_set_user_resuidal_cb(void (*cb)(FLOAT* resuidal, size_t size, void* user_ctx), void* user_ctx, gha_ctx_t ctx);

#endif
