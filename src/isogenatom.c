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

/* Reference (lightest tabulated) isotope mass of an element, or -1. */
static double element_reference_mass(const int atomic_number)
{
    for (int j = 0; j < dlen; j++)
    {
        if (isotope_numbers[j] == atomic_number)
        {
            return isotope_masses[j];
        }
    }
    return -1.0;
}

/* Validate enrichment entries and compute, for one element, the summed
 * heavy fraction and the nominal-mass bin of each heavy isotope. Returns the
 * number of entries for the element, or -1 on invalid input. */
static int collect_enrichment(
    const int atomic_number, const int length,
    const int* enrich_elements, const int* enrich_mass_numbers,
    const double* enrich_fractions, const int enrich_count,
    double* bins, double* fractions, double* total_fraction)
{
    const double reference = element_reference_mass(atomic_number);
    int n = 0;
    *total_fraction = 0.0;
    for (int k = 0; k < enrich_count; k++)
    {
        if (enrich_elements[k] != atomic_number)
        {
            continue;
        }
        const double f = enrich_fractions[k];
        const double heavy_mass =
            atom_isotope_mass(atomic_number, enrich_mass_numbers[k]);
        if (!(f >= 0.0 && f <= 1.0) || heavy_mass < 0.0 || reference < 0.0)
        {
            return -1;
        }
        const double bin = round(heavy_mass - reference);
        if (bin < 0.0 || bin >= length)
        {
            return -1;
        }
        bins[n] = bin;
        fractions[n] = f;
        *total_fraction += f;
        n++;
    }
    if (*total_fraction > 1.0 + 1e-12)
    {
        return -1;
    }
    return n;
}

#define ISOGENATOM_MAX_ENRICH 16

/* FFT of (1 - sum f) * natural + sum f_i * delta(heavy_i). This replaces the
 * precomputed FFT for the element. */
static int setup_enriched_ft(
    const int atomic_number, fftw_complex* outft, const int length,
    const int* enrich_elements, const int* enrich_mass_numbers,
    const double* enrich_fractions, const int enrich_count)
{
    double bins[ISOGENATOM_MAX_ENRICH];
    double fractions[ISOGENATOM_MAX_ENRICH];
    double total = 0.0;
    int n = 0;
    for (int k = 0; k < enrich_count; k++)
    {
        if (enrich_elements[k] == atomic_number) n++;
    }
    if (n > ISOGENATOM_MAX_ENRICH) return -1;
    n = collect_enrichment(
        atomic_number, length, enrich_elements, enrich_mass_numbers,
        enrich_fractions, enrich_count, bins, fractions, &total);
    if (n < 0) return -1;

    double* array = (double*)fftw_malloc((size_t)length * sizeof(double));
    if (array == NULL) return -1;
    memset(array, 0, (size_t)length * sizeof(double));
    construct_isotope_array(atomic_number, array, length);

    double sum = 0.0;
    for (int i = 0; i < length; i++) sum += array[i];
    if (sum <= 0.0)
    {
        fftw_free(array);
        return -1;
    }
    for (int i = 0; i < length; i++) array[i] *= (1.0 - total) / sum;
    for (int k = 0; k < n; k++) array[(int)bins[k]] += fractions[k];

    fftw_plan plan = fftw_plan_dft_r2c_1d(length, array, outft, FFTW_ESTIMATE);
    if (plan == NULL)
    {
        fftw_free(array);
        return -1;
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(array);
    return 0;
}

static int formula_vector_to_probability_dist(
    const int atom_counts[ISOGENATOM_ELEMENT_COUNT],
    const int length,
    float* probability_dist,
    float* max_probability,
    const int* isotope_elements, const int* isotope_mass_numbers,
    const int* isotope_counts, const int isotope_count,
    const int* enrich_elements, const int* enrich_mass_numbers,
    const double* enrich_fractions, const int enrich_count)
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

        int enriched = 0;
        for (int k = 0; k < enrich_count; k++)
        {
            if (enrich_elements[k] == atomic_index + 1) enriched = 1;
        }
        const int setup_status = enriched
            ? setup_enriched_ft(atomic_index + 1, elementft, length,
                                enrich_elements, enrich_mass_numbers,
                                enrich_fractions, enrich_count)
            : setup_ft(atomic_index + 1, elementft, length, ftlen);
        if (setup_status != 0)
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
            &max_probability, NULL, NULL, NULL, 0,
            NULL, NULL, NULL, 0) != 0)
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

float fft_atom_counts_to_dist_enriched(
    const int atom_counts[ISOGENATOM_ELEMENT_COUNT],
    const int* isotope_elements, const int* isotope_mass_numbers,
    const int* isotope_counts, const int isotope_count,
    const int* enrich_elements, const int* enrich_mass_numbers,
    const double* enrich_fractions, const int enrich_count,
    float* isodist, const int isolen, const int offset)
{
    if (atom_counts == NULL || isodist == NULL || isolen <= 0 ||
        offset < 0 || offset >= isolen || isotope_count < 0 || enrich_count < 0 ||
        (isotope_count && (isotope_elements == NULL ||
                          isotope_mass_numbers == NULL || isotope_counts == NULL)) ||
        (enrich_count && (enrich_elements == NULL ||
                         enrich_mass_numbers == NULL || enrich_fractions == NULL))) {
        return -1.0f;
    }
    memset(isodist, 0, (size_t)isolen * sizeof(*isodist));
    for (int k = 0; k < enrich_count; k++) {
        if (enrich_elements[k] < 1 || enrich_elements[k] > ISOGENATOM_ELEMENT_COUNT ||
            atom_counts[enrich_elements[k] - 1] <= 0) {
            return -1.0f;
        }
    }
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
        second /= sum;
        double bins[ISOGENATOM_MAX_ENRICH], fractions[ISOGENATOM_MAX_ENRICH];
        double total = 0.0;
        int n = 0;
        for (int k = 0; k < enrich_count; k++) {
            if (enrich_elements[k] == atomic_index + 1) n++;
        }
        if (n > ISOGENATOM_MAX_ENRICH) return -1.0f;
        n = collect_enrichment(atomic_index + 1, INT_MAX, enrich_elements,
                               enrich_mass_numbers, enrich_fractions, enrich_count,
                               bins, fractions, &total);
        if (n < 0) return -1.0f;
        if (n > 0) {
            first *= 1.0 - total;
            second *= 1.0 - total;
            for (int k = 0; k < n; k++) {
                first += fractions[k] * bins[k];
                second += fractions[k] * bins[k] * bins[k];
            }
        }
        mean += count * first;
        variance += count * fmax(0.0, second - first * first);
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
        isotope_elements, isotope_mass_numbers, isotope_counts, isotope_count,
        enrich_elements, enrich_mass_numbers, enrich_fractions, enrich_count);
    if (status == 0 && max_probability > 0.0f) {
        for (int i = 0; i < isolen - offset; i++) {
            isodist[i + offset] = probabilities[i] / max_probability;
        }
    }
    free(probabilities);
    return status == 0 ? max_probability : -1.0f;
}

float fft_atom_counts_to_dist(
    const int atom_counts[ISOGENATOM_ELEMENT_COUNT],
    const int* isotope_elements, const int* isotope_mass_numbers,
    const int* isotope_counts, const int isotope_count,
    float* isodist, const int isolen, const int offset)
{
    return fft_atom_counts_to_dist_enriched(
        atom_counts, isotope_elements, isotope_mass_numbers, isotope_counts,
        isotope_count, NULL, NULL, NULL, 0, isodist, isolen, offset);
}
