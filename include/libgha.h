#ifndef LIBGHA_H
#define LIBGHA_h

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
 */
gha_ctx_t gha_create_ctx(size_t size);

/*
 * Free GHA context
 */
void gha_free_ctx(gha_ctx_t ctx);

/*
 * Performs one GHA step for given PCM signal,
 * the result will be writen in to given gha_ingo structure
 */
void gha_analyze_one(const FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx);

/*
 * Performs one GHA step and extract analysed harmonic from given PCM signal
 */
void gha_extract_one(FLOAT* pcm, struct gha_info* info, gha_ctx_t ctx);

#endif
