#ifndef LIBGHA_H
#define LIBGHA_h

#define FLOAT float

#include <stddef.h>

typedef struct gha_ctx gha_ctx;

struct gha_info {
	FLOAT freq;
	FLOAT phase;
	FLOAT magnitude;
};

// size must be even
gha_ctx* gha_create_ctx(size_t size);
void gha_free_ctx(gha_ctx* ctx);

// This function performs one GHA step for given PCM signal,
// the result will be writen in to given gha_ingo structure
void gha_analyze_one(const FLOAT* pcm, struct gha_info* info, gha_ctx* ctx);

// Performs one GHA step and extract analysed harmonic from given PCM signal
void gha_extract_one(FLOAT* pcm, struct gha_info* info, gha_ctx* ctx);

#endif
