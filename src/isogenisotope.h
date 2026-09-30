#ifndef ISOGENISOTOPE_H
#define ISOGENISOTOPE_H

#include "isogenatom.h"
#include "fftw3.h"

/* Fixed isotopes are referenced to their own exact mass, so their envelope
 * is a delta at zero. Its FFT is unity at every frequency. This separate
 * precomputed object never changes the natural-abundance fftarray. */
extern const double fixed_isotope_fft[33][2];
ISOGENATOM_EXPORTS double atom_isotope_mass(int atomic_number, int mass_number);
int setup_fixed_isotope_ft(int atomic_number, int mass_number,
                          fftw_complex* output, int length);

#endif
