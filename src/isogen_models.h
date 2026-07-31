#ifndef ISOGEN_MODELS_H
#define ISOGEN_MODELS_H

#include "isogendep.h"

/*
 * Load weights from the raw binary format produced by
 * save_model_to_binary(). The caller must initialize weights with
 * SetupWeights() using dimensions that match the file.
 *
 * Returns 0 on success and -1 if the arguments are invalid, the file cannot
 * be read, or its size does not exactly match the allocated model. On failure,
 * weights are left unchanged and remain owned by the caller.
 */
ISOGENDEP_EXPORTS int LoadWeightsFromFile(struct IsoGenWeights weights, const char *filename);

/* Shared implementation used by the public custom-model prediction APIs. */
EXTERN float isogen_model_to_dist_from_file(const float *vector, int vector_length,
                                            float *isodist, int isolen, int offset,
                                            const char *filename);

#endif /* ISOGEN_MODELS_H */
