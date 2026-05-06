#include <fstream>
#include "isogenpep.h"


void print_isodist(const float *isodist, const int size) {
    // Make a command line printout of the spectrum
    for (int i=10; i>=0; i--) {
        for (int j=0; j<size; j++) {
            const float val = isodist[j];
            if ((float) i < val*10) {
                printf("*");
            }
            else {
                printf(" ");
            }
        }
        printf("\n");
    }

    // Print floats
    // for (int i=0; i<size; i++) {
    //     printf("%.2f ", isodist[i]);
    // }
    printf("\n");
}


int main(const int argc, char *argv[])
{
    printf("Running IsoGen Exe\n");
    // Get arguments from command line
    if (argc < 2)
    {
        printf("Usage: %s -mass <mass float>\n", argv[0]);
        return 1;
    }


    if (argc == 2) {
        printf("Testing IsoGen from exe\n");
        char* filename = argv[1];
        printf("Testing IsoGen on %s\n", filename);

        run_file(filename);
    }

    if (argc == 3) {
        // Check if argv[1] is -mass
        if (strcmp(argv[1], "-mass") == 0) {
            // Run isogenmass
            float mass = 0.0f;
            mass = strtof(argv[2], nullptr);
            int isolen = nn_pep_mass_to_isolen(mass);
            float *isodist = (float *) calloc(isolen, sizeof(float));
            if(isodist == nullptr) {
                printf("Error allocating memory for isodist\n");
                return 1;
            }
            printf("Running NN Peptide Mass Prediction on: %.2f Da\n", mass);
            nn_pep_mass_to_dist(mass, isodist, isolen, 0);
            print_isodist(isodist, isolen);
            free(isodist);
        }
    }
    return 0;

}
