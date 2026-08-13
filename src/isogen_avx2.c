#include <immintrin.h>

#include "isogendep.h"


void matrix_vector_multiply_avx2(const float* matrix, const float* vector, const float* bias,
                                 float* result, const int N, const int M,
                                 const int ss_flag, const int sigmoid_flag)
{
    for (int i = 0; i < N; i++)
    {
        float value = bias[i];
        __m256 sum = _mm256_setzero_ps();
        int j = 0;
        for (; j <= M - 8; j += 8)
        {
            const __m256 matrix_values = _mm256_loadu_ps(&matrix[i * M + j]);
            const __m256 vector_values = _mm256_loadu_ps(&vector[j]);
            sum = _mm256_fmadd_ps(matrix_values, vector_values, sum);
        }

        float partial_sums[8];
        _mm256_storeu_ps(partial_sums, sum);
        value += partial_sums[0] + partial_sums[1] + partial_sums[2] + partial_sums[3] +
                 partial_sums[4] + partial_sums[5] + partial_sums[6] + partial_sums[7];

        for (; j < M; j++)
        {
            value += matrix[i * M + j] * vector[j];
        }
        if (ss_flag)
        {
            value = softsign(value);
        }
        if (sigmoid_flag)
        {
            value = sigmoid(value);
        }
        result[i] = value;
    }
}
