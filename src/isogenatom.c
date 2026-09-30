#include "isogenatom.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fftw3.h"
#include "isogendep.h"
#include "isogenisotope.h"

extern const int dlen;

static int symbol_to_atomic_number(const char* symbol, const size_t symbol_length)
{
    for (int i = 0; i < ISOGENATOM_ELEMENT_COUNT; i++)
    {
        if (strlen(element_names[i]) == symbol_length &&
            strncmp(element_names[i], symbol, symbol_length) == 0)
        {
            return i + 1;
        }
    }
    return -1;
}

int atom_formula_to_vector(
    const char* formula,
    int atom_counts[ISOGENATOM_ELEMENT_COUNT])
{
    if (formula == NULL || atom_counts == NULL)
    {
        return -1;
    }

    memset(atom_counts, 0, ISOGENATOM_ELEMENT_COUNT * sizeof(*atom_counts));

    const char* cursor = formula;
    int total_atoms = 0;
    int found_element = 0;

    while (*cursor != '\0')
    {
        if (isspace((unsigned char)*cursor))
        {
            cursor++;
            continue;
        }

        if (!isupper((unsigned char)*cursor))
        {
            return -1;
        }

        const char* symbol = cursor++;
        size_t symbol_length = 1;
        if (islower((unsigned char)*cursor))
        {
            cursor++;
            symbol_length++;
        }

        const int atomic_number = symbol_to_atomic_number(symbol, symbol_length);
        if (atomic_number < 1)
        {
            return -1;
        }

        int count = 0;
        int has_count = 0;
        while (isdigit((unsigned char)*cursor))
        {
            const int digit = *cursor - '0';
            if (count > (INT_MAX - digit) / 10)
            {
                return -1;
            }
            count = count * 10 + digit;
            has_count = 1;
            cursor++;
        }
        if (!has_count)
        {
            count = 1;
        }

        const int index = atomic_number - 1;
        if (atom_counts[index] > INT_MAX - count ||
            total_atoms > INT_MAX - count)
        {
            return -1;
        }
        atom_counts[index] += count;
        total_atoms += count;
        found_element = 1;
    }

    return found_element && total_atoms > 0 ? 0 : -1;
}

static int formula_vector_to_probability_dist(
    const int atom_counts[ISOGENATOM_ELEMENT_COUNT],
    const int length,
    float* probability_dist,
    float* max_probability,
    const int* isotope_elements, const int* isotope_mass_numbers,
    const int* isotope_counts, const int isotope_count)
{
    const int ftlen = length / 2 + 1;
    fftw_complex* allft =
        (fftw_complex*)fftw_malloc((size_t)ftlen * sizeof(fftw_complex));
    fftw_complex* elementft =
        (fftw_complex*)fftw_malloc((size_t)ftlen * sizeof(fftw_complex));
    double* buffer = (double*)fftw_malloc((size_t)length * sizeof(double));

    if (allft == NULL || elementft == NULL || buffer == NULL)
    {
        fftw_free(allft);
        fftw_free(elementft);
        fftw_free(buffer);
        return -1;
    }

    for (int i = 0; i < ftlen; i++)
    {
        allft[i][0] = 1.0;
        allft[i][1] = 0.0;
    }

    for (int atomic_index = 0;
         atomic_index < ISOGENATOM_ELEMENT_COUNT;
         atomic_index++)
    {
        const int count = atom_counts[atomic_index];
        if (count == 0)
        {
            continue;
        }

        if (setup_ft(atomic_index + 1, elementft, length, ftlen) != 0)
        {
            fftw_free(allft);
            fftw_free(elementft);
            fftw_free(buffer);
            return -1;
        }

        for (int i = 0; i < ftlen; i++)
        {
            double powered_real;
            double powered_imag;
            double product_real;
            double product_imag;

            complex_power(
                elementft[i],
                count,
                &powered_real,
                &powered_imag
            );
            complex_multiplication(
                allft[i][0],
                allft[i][1],
                powered_real,
                powered_imag,
                &product_real,
                &product_imag
            );
            allft[i][0] = product_real;
            allft[i][1] = product_imag;
        }
    }

    for (int j = 0; j < isotope_count; j++) {
        if (isotope_counts[j] < 0 || setup_fixed_isotope_ft(
                isotope_elements[j], isotope_mass_numbers[j], elementft, length) != 0) {
            fftw_free(allft);
            fftw_free(elementft);
            fftw_free(buffer);
            return -1;
        }
        /* The fixed-isotope FFT is unity relative to the labeled mass origin.
         * Explicitly combine it without altering the natural isotope tables. */
        for (int i = 0; i < ftlen; i++) {
            allft[i][0] *= elementft[i][0];
            allft[i][1] *= elementft[i][0];
        }
    }

    fftw_plan inverse_plan =
        fftw_plan_dft_c2r_1d(length, allft, buffer, FFTW_ESTIMATE);
    if (inverse_plan == NULL)
    {
        fftw_free(allft);
        fftw_free(elementft);
        fftw_free(buffer);
        return -1;
    }

    fftw_execute(inverse_plan);
    fftw_destroy_plan(inverse_plan);
    fftw_free(allft);
    fftw_free(elementft);

    *max_probability = (float)normalize_isodist(buffer, length);
    for (int i = 0; i < length; i++)
    {
        probability_dist[i] = (float)buffer[i];
    }
    fftw_free(buffer);
    return 0;
}

float fft_atom_formula_to_dist(
    const char* formula,
    float* isodist,
    const int isolen,
    const int offset)
{
    if (isodist == NULL || isolen <= 0 || offset < 0 || offset >= isolen)
    {
        return -1.0f;
    }

    for (int i = 0; i < isolen; i++)
    {
        isodist[i] = 0.0f;
    }

    int atom_counts[ISOGENATOM_ELEMENT_COUNT];
    if (atom_formula_to_vector(formula, atom_counts) != 0)
    {
        return -1.0f;
    }

    float* probability_dist =
        (float*)calloc((size_t)isolen, sizeof(float));
    if (probability_dist == NULL)
    {
        return -1.0f;
    }

    float max_probability = 0.0f;
    if (formula_vector_to_probability_dist(
            atom_counts,
            isolen,
            probability_dist,
            &max_probability, NULL, NULL, NULL, 0) != 0)
    {
        free(probability_dist);
        return -1.0f;
    }

    if (max_probability > 0.0f)
    {
        const int copy_length = isolen - offset;
        for (int i = 0; i < copy_length; i++)
        {
            isodist[i + offset] = probability_dist[i] / max_probability;
        }
    }

    free(probability_dist);
    return max_probability;
}

void isogen_atom(const char* formula, float* isodist, const int isolen)
{
    (void)fft_atom_formula_to_dist(formula, isodist, isolen, 0);
}

float fft_atom_counts_to_dist(
    const int atom_counts[ISOGENATOM_ELEMENT_COUNT],
    const int* isotope_elements, const int* isotope_mass_numbers,
    const int* isotope_counts, const int isotope_count,
    float* isodist, const int isolen, const int offset)
{
    if (atom_counts == NULL || isodist == NULL || isolen <= 0 ||
        offset < 0 || offset >= isolen || isotope_count < 0 ||
        (isotope_count && (isotope_elements == NULL ||
                          isotope_mass_numbers == NULL || isotope_counts == NULL))) {
        return -1.0f;
    }
    memset(isodist, 0, (size_t)isolen * sizeof(*isodist));
    double mean = 0.0;
    double variance = 0.0;
    for (int atomic_index = 0; atomic_index < ISOGENATOM_ELEMENT_COUNT; atomic_index++) {
        const int count = atom_counts[atomic_index];
        if (count < 0) return -1.0f;
        if (!count) continue;
        double reference = -1.0;
        double sum = 0.0, first = 0.0, second = 0.0;
        for (int j = 0; j < dlen; j++) {
            if (isotope_numbers[j] != atomic_index + 1) continue;
            if (reference < 0) reference = isotope_masses[j];
            const double shift = round(isotope_masses[j] - reference);
            const double p = isotope_abundances[j];
            sum += p;
            first += p * shift;
            second += p * shift * shift;
        }
        if (sum <= 0) return -1.0f;
        first /= sum;
        mean += count * first;
        variance += count * fmax(0.0, second / sum - first * first);
    }
    /* Zero-pad before convolution: short requested outputs must truncate,
     * rather than fold the envelope tail into the first isotope bins. */
    const double needed = fmax((double)isolen, mean + 12.0 * sqrt(variance) + 32.0);
    int length = 64;
    while (length < needed) {
        if (length > INT_MAX / 2) return -1.0f;
        length *= 2;
    }
    float* probabilities = (float*)calloc((size_t)length, sizeof(float));
    if (probabilities == NULL) return -1.0f;
    float max_probability = 0.0f;
    const int status = formula_vector_to_probability_dist(
        atom_counts, length, probabilities, &max_probability,
        isotope_elements, isotope_mass_numbers, isotope_counts, isotope_count);
    if (status == 0 && max_probability > 0.0f) {
        for (int i = 0; i < isolen - offset; i++) {
            isodist[i + offset] = probabilities[i] / max_probability;
        }
    }
    free(probabilities);
    return status == 0 ? max_probability : -1.0f;
}
