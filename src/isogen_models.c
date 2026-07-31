#include "isogen_models.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int weights_have_valid_storage(const struct IsoGenWeights weights) {
    if (weights.tot <= 0 || weights.ml1 <= 0 || weights.vl2 <= 0 ||
        weights.ml2 <= 0 || weights.vl3 <= 0 || weights.ml3 <= 0 ||
        weights.vl4 <= 0 || weights.w1 == NULL || weights.b1 == NULL ||
        weights.w2 == NULL || weights.b2 == NULL ||
        weights.w3 == NULL || weights.b3 == NULL) {
        return 0;
    }

    const int lengths[] = {
        weights.ml1, weights.vl2, weights.ml2,
        weights.vl3, weights.ml3, weights.vl4
    };
    size_t total = 0;
    for (size_t index = 0; index < sizeof(lengths) / sizeof(lengths[0]); ++index) {
        if ((size_t)lengths[index] > (size_t)weights.tot - total) {
            return 0;
        }
        total += (size_t)lengths[index];
    }
    return total == (size_t)weights.tot;
}

int LoadWeightsFromFile(const struct IsoGenWeights weights, const char *filename) {
    if (filename == NULL || filename[0] == '\0' || !weights_have_valid_storage(weights)) {
        return -1;
    }

    const size_t expected_size = (size_t)weights.tot * sizeof(float);
    unsigned char *model_weights = (unsigned char *)malloc(expected_size);
    if (model_weights == NULL) {
        return -1;
    }

    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        free(model_weights);
        return -1;
    }

    const size_t bytes_read = fread(model_weights, 1, expected_size, file);
    unsigned char extra_byte;
    const size_t extra_bytes_read = fread(&extra_byte, 1, 1, file);
    const int read_failed = ferror(file);
    fclose(file);

    if (bytes_read != expected_size || extra_bytes_read != 0 || read_failed) {
        free(model_weights);
        return -1;
    }

    LoadWeights(weights, model_weights);
    free(model_weights);
    return 0;
}

float isogen_model_to_dist_from_file(const float *vector, const int vector_length,
                                     float *isodist, const int isolen, const int offset,
                                     const char *filename) {
    if (vector == NULL || isodist == NULL || vector_length <= 0 || isolen <= 0 ||
        offset < 0 || offset >= isolen || filename == NULL || filename[0] == '\0') {
        return -1.0f;
    }

    memset(isodist, 0, (size_t)isolen * sizeof(*isodist));

    struct IsoGenWeights weights = SetupWeights(vector_length, isolen);
    if (!weights_have_valid_storage(weights)) {
        FreeIsogenWeights(weights);
        return -1.0f;
    }

    if (LoadWeightsFromFile(weights, filename) != 0) {
        FreeIsogenWeights(weights);
        return -1.0f;
    }

    float *nn_isodist = (float *)calloc((size_t)isolen, sizeof(*nn_isodist));
    if (nn_isodist == NULL) {
        FreeIsogenWeights(weights);
        return -1.0f;
    }

    if (neural_net(vector, nn_isodist, weights) != 0) {
        free(nn_isodist);
        FreeIsogenWeights(weights);
        return -1.0f;
    }
    FreeIsogenWeights(weights);

    for (int index = 0; index < isolen - offset; ++index) {
        isodist[index + offset] = nn_isodist[index];
    }
    free(nn_isodist);

    float maxval = 0.0f;
    for (int index = 0; index < isolen; ++index) {
        if (isodist[index] > maxval) {
            maxval = isodist[index];
        }
    }
    if (maxval > 0.0f) {
        for (int index = 0; index < isolen; ++index) {
            isodist[index] /= maxval;
        }
    }
    return maxval;
}
