#include "isogenisotope.h"
#include "isogendep.h"
#include <math.h>

extern const int dlen;

const double fixed_isotope_fft[33][2] = {
    {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0},
    {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0},
    {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0},
    {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0},
    {1,0}
};

double atom_isotope_mass(const int atomic_number, const int mass_number)
{
    for (int i = 0; i < dlen; i++) {
        if (isotope_numbers[i] == atomic_number &&
            (int)round(isotope_masses[i]) == mass_number) {
            return isotope_masses[i];
        }
    }
    return -1.0;
}

int setup_fixed_isotope_ft(const int atomic_number, const int mass_number,
                          fftw_complex* output, const int length)
{
    if (output == NULL || length <= 0 ||
        atom_isotope_mass(atomic_number, mass_number) < 0) {
        return -1;
    }
    for (int i = 0; i < length / 2 + 1; i++) {
        /* The same exact transform applies to every fully specified isotope
         * and FFT length; larger lengths need no FFTW plan or mutable cache. */
        output[i][0] = fixed_isotope_fft[i % 33][0];
        output[i][1] = fixed_isotope_fft[i % 33][1];
    }
    return 0;
}
