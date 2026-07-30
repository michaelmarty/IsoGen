#ifndef ISOGENATOM_H
#define ISOGENATOM_H

#ifdef __cplusplus
    #define ISOGENATOM_EXTERN extern "C"
#else
    #define ISOGENATOM_EXTERN
#endif

#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENATOM_EXPORTS ISOGENATOM_EXTERN __declspec(dllexport)
    #else
        #define ISOGENATOM_EXPORTS ISOGENATOM_EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENATOM_EXPORTS ISOGENATOM_EXTERN __declspec(dllimport)
    #else
        #define ISOGENATOM_EXPORTS ISOGENATOM_EXTERN
    #endif
#endif

#define ISOGENATOM_ELEMENT_COUNT 109

/*
 * Convert a formula such as "C6H12O6" into counts indexed by atomic
 * number minus one. Returns 0 on success and -1 for an invalid formula.
 */
ISOGENATOM_EXPORTS int atom_formula_to_vector(
    const char* formula,
    int atom_counts[ISOGENATOM_ELEMENT_COUNT]
);

/*
 * Calculate a base-peak-normalized isotope distribution for an elemental
 * formula. The returned value is the maximum probability before base-peak
 * normalization, matching the peptide and RNA FFT interfaces. A negative
 * value indicates invalid input or an allocation failure.
 */
ISOGENATOM_EXPORTS float fft_atom_formula_to_dist(
    const char* formula,
    float* isodist,
    int isolen,
    int offset
);

/* Compatibility entry point retained from the original public declaration. */
ISOGENATOM_EXPORTS void isogen_atom(
    const char* formula,
    float* isodist,
    int isolen
);

#endif /* ISOGENATOM_H */
