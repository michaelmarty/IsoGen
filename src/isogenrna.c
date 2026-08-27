#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogen_models.h"
#include "isogen_rnaveragine_model32.h"
#include "isogen_rnaveragine_model64.h"
#include "isogen_rnaveragine_model128.h"
#include "isogenrna_model_64.h"
#include "isogenrna_model_128.h"



const char rnaOrder[] = "ACGU";

//All rna averagine values assume 1 phospho per nucleotide
const float rnaveragine_mass = 320.283814;
const double rnaveragine_comp_numerical[] = {9.50, 10.75, 3.75, 7.0, 0};

static const unsigned char nt_lookup[256] = {
    ['A']=1, ['C']=2, ['G']=3, ['U']=4,
    ['a']=1, ['c']=2, ['g']=3, ['u']=4
};

const int rna_vectors[][5] = {
    {10,11,5,6,0},
    {9,11,3,7,0},
    {10,11,5,7,0},
    {9,10,2,8,0}
};

const double rnaveragine_coeff[] = {0.029661, 0.033564, 0.01171, 0.02186, 0};

static int prepare_isodist_output(float* isodist, const int isolen, const int offset)
{
    if (isodist == NULL || isolen <= 0 || offset < 0 || offset >= isolen)
    {
        return -1;
    }
    memset(isodist, 0, (size_t)isolen * sizeof(*isodist));
    return 0;
}

static int nn_output_ready(const struct IsoGenWeights weights, const float* nn_isodist)
{
    return nn_isodist != NULL &&
           weights.w1 != NULL && weights.b1 != NULL &&
           weights.w2 != NULL && weights.b2 != NULL &&
           weights.w3 != NULL && weights.b3 != NULL;
}

int fft_rna_len_to_isolen(const int rna_len);
int fft_rna_mass_to_isolen(float mass);


void rna_mass_to_list(float initialMass, int* fftlist)
//Mass -> List 5 length list of number of {C, H, N, O, S}
{
    float rnaMass = rnaveragine_mass;
    //Addition of an extra phospho here simplifies the monomer number calculation
    float valuePer = (initialMass + 95.9534) / rnaMass;

    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(rnaveragine_comp_numerical[i] * valuePer);
    }

    //Correct for additional one less phospho than monomers
    fftlist[3] -= 4;
}



float fft_rna_mass_to_dist(float mass, float* isodist, int isolen, int offset)
// Mass to dist
{
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    int* fftlist = (int*)calloc(5, sizeof(int));
    if (fftlist == NULL) {
        return -1.0f;
    }
    rna_mass_to_list(mass, fftlist);

    int fft_isolen = fft_rna_mass_to_isolen(mass);
    float* fft_isodist = (float*)calloc(fft_isolen, sizeof(float));
    if (fft_isodist == NULL) {
        free(fftlist);
        return -1.0f;
    }

    float max_val = fft_list_to_dist(fftlist, fft_isolen, fft_isodist);
    if (max_val < 0.0f) {
        free(fft_isodist);
        free(fftlist);
        return -1.0f;
    }

    int copy_len = fft_isolen < isolen - offset ? fft_isolen : isolen - offset;
    for (int i = copy_len - 1; i >= 0; i--)
    {
        isodist[i + offset] = fft_isodist[i];
    }

    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }
    free(fft_isodist);
    free(fftlist);
    return max_val;
}


int nt_to_index(char nt) {
    for (int i = 0;i<4;i++) {
        if (rnaOrder[i] == nt) {
            return i;
        }
    }
    return -1;
}



int fft_rna_len_to_isolen(const int rna_len) {
    if (rna_len < 200) {
        return 64;
    }
    if (rna_len < 501) {
        return 128;
    }
    return 1024;
}

int nn_rna_mass_to_isolen(const float mass) {
    if (mass < 24000) {
        return 32;
    }
    if (mass < 65000) {
        return 64;
    }
    if (mass < 165000) {
        return 128;
    }
    return -1;
}

int fft_rna_mass_to_isolen(float mass) {
    if (mass < 65000) {
        return 64;
    }
    if (mass < 165000) {
        return 128;
    }
    return 1024;
}


int rna_seq_to_vector(const char* seq, float* vector)
//RNA -> Vector 4 length list of number of A C G U
{
    int len = strlen(seq);
    for (int i = 0; i < len; i++)
    {
        int nt_index = nt_to_index(seq[i]);
        if (nt_index == -1) {
            printf("Unexpected nucleotide in sequence: %c\n", seq[i]);
        }
        else {
            vector[nt_index] += 1.0f;
        }
    }
    return len;
}


int rna_seq_to_fftlist(const char* sequence, int* fftlist)
{
    if (sequence == NULL || fftlist == NULL) {
        return -1;
    }
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 0; // Carbon
    fftlist[1] = 0; // Hydrogen
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 0; // Oxygen
    fftlist[4] = 0; // Sulfur


    int seq_len = strlen(sequence);
    for (int i = 0; i < seq_len; i++) {
        int nt_index = nt_to_index(sequence[i]);
        if (nt_index == -1) {
            printf("Unexpected nucleotide in sequence:%c\n", sequence[i]);
        }
        else {
            for (int j = 0;j<5;j++) {
                fftlist[j] += rna_vectors[nt_index][j];
            }
        }
    }

    return seq_len;
}


//fft
float fft_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    int* fftlist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (fftlist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return -1.0f;
    }
    int rna_len = rna_seq_to_fftlist(sequence, fftlist);
    if (rna_len < 0) {
        free(fftlist);
        return -1.0f;
    }

    int fft_isolen = fft_rna_len_to_isolen(rna_len);
    float* fft_isodist = (float*)calloc(fft_isolen, sizeof(float));
    if (fft_isodist == NULL) {
        free(fftlist);
        return -1.0f;
    }

    float maxval = fft_list_to_dist(fftlist, fft_isolen, fft_isodist);
    free(fftlist);
    if (maxval < 0.0f) {
        free(fft_isodist);
        return -1.0f;
    }

    int copy_len = fft_isolen < isolen - offset ? fft_isolen : isolen - offset;
    for (int i = copy_len - 1; i >= 0; i--)
    {
        isodist[i + offset] = fft_isodist[i];
    }

    if (maxval > 0.0f) {
        for (int i = 0;i<isolen;i++) {
            isodist[i] /= maxval;
        }
    }

    free(fft_isodist);
    return maxval;
}


float nn_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        return -1.0f;
    }

    mass_to_vector(mass, vector);

    int nn_isolen = nn_rna_mass_to_isolen(mass);

    if (nn_isolen == -1) {
        printf("Error: Mass outside of allowed NN mass range: %f\n", mass);
        free(vector);
        return -1.0f;
    }

    float* nn_isodist = (float*)calloc(nn_isolen, sizeof(float));

    struct IsoGenWeights weights = SetupWeights(5, nn_isolen);
    if (nn_isolen == 32){ weights = LoadWeights(weights, isogen_rnaveragine_model32_bin); }
    else if ( nn_isolen == 64 ){ weights = LoadWeights(weights, isogen_rnaveragine_model64_bin); }
    else { weights = LoadWeights(weights, isogen_rnaveragine_model128_bin); }

    if (!nn_output_ready(weights, nn_isodist)) {
        free(vector);
        free(nn_isodist);
        FreeIsogenWeights(weights);
        return -1.0f;
    }

    if (neural_net(vector, nn_isodist, weights) != 0) {
        free(vector);
        free(nn_isodist);
        FreeIsogenWeights(weights);
        return -1.0f;
    }
    free(vector);
    FreeIsogenWeights(weights);

    if (nn_isolen < isolen) {
        for (int i = nn_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    free(nn_isodist);

    float maxval = 0.0f;
    for (int i = 0; i < isolen; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }
    return maxval;
}

float nn_rna_mass_to_dist_custom(const float mass, float* isodist, const int isolen, const int offset,
                                 const char* model_path) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }

    float vector[5] = {0};
    mass_to_vector(mass, vector);
    return isogen_model_to_dist_from_file(vector, 5, isodist, isolen, offset, model_path);
}


float nn_rna_seq_to_dist(const char* seq, float* isodist, int isolen, int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    float* vector = calloc(4, sizeof(float));
    if (vector == NULL) {
        return -1.0f;
    }

    int len = rna_seq_to_vector(seq, vector);

    struct IsoGenWeights weights = {0};
    float* nn_isodist = NULL;
    int nn_isolen = 0;

    if (len >= 1 && len <= 200) {
        weights = SetupWeights(4, 64);
        weights = LoadWeights(weights, isogenrna_model_64_bin);
        nn_isolen = 64;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    if (len >= 201 && len <= 500) {
        weights = SetupWeights(4, 128);
        weights = LoadWeights(weights, isogenrna_model_128_bin);
        nn_isolen = 128;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    if (len > 500) {
        printf("Error: Sequence length outside of allowed NN range: %i\n", len);
        free(vector);
        return -1.0f;
    }
    if (len < 1) {
        free(vector);
        return -1.0f;
    }

    if (!nn_output_ready(weights, nn_isodist)) {
        free(vector);
        free(nn_isodist);
        FreeIsogenWeights(weights);
        return -1.0f;
    }

    if (neural_net(vector, nn_isodist, weights) != 0) {
        free(vector);
        free(nn_isodist);
        FreeIsogenWeights(weights);
        return -1.0f;
    }
    free(vector);
    FreeIsogenWeights(weights);

    if (nn_isolen < isolen) {
        for (int i = nn_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    free(nn_isodist);


    float maxval = 0.0f;
    for (int i = 0; i < isolen; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }

    if (maxval > 0.0f) {
        for (int i = 0;i< isolen; i++) {
            isodist[i] /= maxval;
        }
    }
    return maxval;
}

float nn_rna_seq_to_dist_custom(const char* seq, float* isodist, const int isolen, const int offset,
                                const char* model_path) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0 || seq == NULL || seq[0] == '\0') {
        return -1.0f;
    }

    float vector[4] = {0};
    rna_seq_to_vector(seq, vector);
    return isogen_model_to_dist_from_file(vector, 4, isodist, isolen, offset, model_path);
}

//BRAIN
float brain_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0 || !isfinite(mass) || mass <= 0.0f) {
        return -1.0f;
    }
    int* brain_list = (int*)calloc(5, sizeof(int));
    if (brain_list == NULL) {
        return -1.0f;
    }
    rna_mass_to_list(mass, brain_list);

    int brain_isolen = fft_rna_mass_to_isolen(mass);

    float* brain_isodist = (float*)calloc(brain_isolen, sizeof(float));
    if (brain_isodist == NULL) {
        free(brain_list);
        return -1.0f;
    }

    float max_val = brain_list_to_dist(brain_list, brain_isolen, brain_isodist);
    if (max_val < 0.0f || !isfinite(max_val)) {
        free(brain_isodist);
        free(brain_list);
        return -1.0f;
    }

    if (brain_isolen < isolen) {
        for (int i = brain_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(brain_isodist);
    free(brain_list);

    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }

    return max_val;
}

//BRAIN
float brain_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0 || sequence == NULL || sequence[0] == '\0') {
        return -1.0f;
    }
    int* formulalist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return -1.0f;
    }
    int len = rna_seq_to_fftlist(sequence, formulalist);
    if (len < 1) {
        free(formulalist);
        return -1.0f;
    }

    int brain_isolen = fft_rna_len_to_isolen(len);

    float* brain_isodist = (float*)calloc(brain_isolen, sizeof(float));
    if (brain_isodist == NULL) {
        free(formulalist);
        return -1.0f;
    }

    float maxval = brain_list_to_dist(formulalist, brain_isolen, brain_isodist);
    if (maxval < 0.0f || !isfinite(maxval)) {
        free(brain_isodist);
        free(formulalist);
        return -1.0f;
    }

    if (brain_isolen < isolen) {
        for (int i = brain_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(brain_isodist);
    free(formulalist);

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }

    return maxval;
}
