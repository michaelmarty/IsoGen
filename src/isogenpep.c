#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "isogendep.h"
#include "isogen_models.h"
#include "isogenpep.h"
#include "isogenpep_model_16.h"
#include "isogenpep_model_64.h"
#include "isogenpep_model_128.h"
#include "isogenmass_model_8.h"
#include "isogenmass_model_32.h"
#include "isogenmass_model_64.h"
#include "isogenmass_model_128.h"


double massavgine = 111.1254;
double avgine[5] = {4.9384, 7.7583, 1.3577, 1.4773, 0.0417};
double averagine_coeffs[5] = {0.044440, 0.069816, 0.012218, 0.013294, 0.0003753};
int numaminoacids = 23;

int num_simp_elements = 5;

//{C, H, N, O, S}
const int aa_vectors[][5] = {
    {3,5,1,1,0},
    {3,5,1,1,1},
    {4,5,1,3,0},
    {5,7,2,3,0},
    {9,9,1,1,0},
    {2,3,1,1,0},
    {6,7,3,1,0},
    {6,11,1,1,0},
    {6,12,2,1,0},
    {6,11,1,1,0},
    {5,9,1,1,1},
    {4,3,2,2,0},
    {5,7,1,1,0},
    {5,7,2,2,0},
    {6,12,4,1,0},
    {5,7,1,2,0},
    {4,7,1,2,0},
    {5,9,1,1,0},
    {11,10,2,1,0},
    {9,9,1,2,0}
};

static const unsigned char aa_lookup[256] = {
    ['A']=1, ['C']=2, ['D']=3, ['E']=4, ['F']=5, ['G']=6, ['H']=7, ['I']=8, ['K']=9, ['L']=10, ['M']=11,
    ['N']=12, ['P']=13, ['Q']=14, ['R']=15, ['S']=16, ['T']=17, ['V']=18, ['W']=19, ['Y']=20,
    ['a']=1, ['c']=2, ['d']=3, ['e']=4, ['f']=5, ['g']=6, ['h']=7, ['i']=8, ['k']=9, ['l']=10, ['m']=11,
    ['n']=12, ['p']=13, ['q']=14, ['r']=15, ['s']=16, ['t']=17, ['v']=18, ['w']=19, ['y']=20
};

const char aaorder[] = "ACDEFGHIKLMNPQRSTVWY";

const char *pep_encoding_elements[] = {"C", "H", "N", "O", "S"};

#define ISO_LEN 32

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

// Isotope Parameters
float isoparams[10] = {
    1.00840852e+00f, 1.25318718e-03f, 2.37226341e+00f, 8.19178000e-04f, -4.37741951e-01f,
    6.64992972e-04f, 9.94230511e-01f, 4.64975237e-01f, 1.00529041e-02f, 5.81240792e-01f
};

// Function used in isotope distribution calculation
float isotopemid(const float mass, const float *isoparams) {
    const float a = isoparams[4];
    const float b = isoparams[5];
    const float c = isoparams[6];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopesig(const float mass, const float *isoparams) {
    const float a = isoparams[7];
    const float b = isoparams[8];
    const float c = isoparams[9];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopealpha(const float mass, const float *isoparams) {
    const float a = isoparams[0];
    const float b = isoparams[1];
    return a * expf(-mass * b);
}

// Function used in isotope distribution calculation
float isotopebeta(const float mass, const float *isoparams) {
    const float a = isoparams[2];
    const float b = isoparams[3];
    return a * expf(-mass * b);
}

// Averagine distribution generation by curve fitting.
float pep_mass_to_dist_fitting(const float mass, float * isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    const float mid = isotopemid(mass, isoparams);
    const float sig = isotopesig(mass, isoparams);
    if (sig <= 0.0f) {
        printf("Error: Sigma isotope parameter must be positive\n");
        return -1.0f;
    }
    const float alpha = isotopealpha(mass, isoparams);
    const float amp = 1.0f - alpha;
    const float beta = isotopebeta(mass, isoparams);

    float max = 0.0f;
    for (int k = 0; k < isolen-offset; k++) {
        const float kfloat = (float) k;
        const float e = alpha * expf(-kfloat * beta);
        const float g = amp / (sig * 2.50662827f) * expf(-powf(kfloat - mid, 2) / (2 * powf(sig, 2)));
        const float val = e + g;
        if(val> max) {max = val;}
        isodist[k+offset] = val;
    }
    return max;
}


int aa_to_index(char aa) {
    for (int i = 0;i<20;i++) {
        if (aa == aaorder[i]){return i;}
    }
    return -1;
}

int pep_seq_to_aacount(const char* seq) {
    int length = strlen(seq);

    int aas = 0;

    int in_mod = 0;

    for (int i = 0; i < length; i++) {
        if (in_mod == 0) {
            if (seq[i] == '[') {
                in_mod = 1;
                continue;
            }
            aas += 1;
        }
        else {
            if (seq[i] == ']'){ in_mod = 0; }
        }
    }
    return aas;
}


//NN
int pep_seq_to_nnvector(const char* seq, float* vector) {
    int length = strlen(seq);

    int aas = 0;

    int in_mod = 0;

    for (int i = 0; i < length; i++) {
        if (in_mod == 0) {
            if (seq[i] == '[') {
                in_mod = 1;
                continue;
            }
            aas += 1;
            int aaindex = aa_to_index(seq[i]);
            if (aaindex != -1){ vector[aaindex] += 1.0f; }
        }
        else {
            if (seq[i] == ']'){ in_mod = 0; }
        }
    }
    return aas;
}


float nn_pep_seq_to_dist(const char* seq, float* isodist, int isolen, int offset){
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    float* vector = (float *) calloc(20, sizeof(float));
    if (vector == NULL) {
        return -1.0f;
    }
    int aas = pep_seq_to_nnvector(seq, vector);

    if (aas > 1000) {
        printf("Sequence contains too many amino acids (>1000).");
        free(vector);
        return -1.0f;
    }

    struct IsoGenWeights weights = {0};
    float* nn_isodist = NULL;
    int nn_isolen = 0;

    if (aas >= 1 && aas <= 50) {
        weights = SetupWeights(20, 16);
        weights = LoadWeights(weights, isogenpep_model_16_bin);
        nn_isolen = 16;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    else if (aas > 50 && aas <= 300) {
        weights = SetupWeights(20, 64);
        weights = LoadWeights(weights, isogenpep_model_64_bin);
        nn_isolen = 64;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    else {
        weights = SetupWeights(20, 128);
        weights = LoadWeights(weights, isogenpep_model_128_bin);
        nn_isolen = 128;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
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

float nn_pep_seq_to_dist_custom(const char* seq, float* isodist, const int isolen, const int offset,
                                const char* model_path) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0 || seq == NULL || seq[0] == '\0') {
        return -1.0f;
    }

    float *vector = (float *)calloc(20, sizeof(*vector));
    if (vector == NULL) {
        return -1.0f;
    }
    pep_seq_to_nnvector(seq, vector);

    const float maxval = isogen_model_to_dist_from_file(
        vector, 20, isodist, isolen, offset, model_path);
    free(vector);
    return maxval;
}

//fft
void add_mod_to_fftlist(const char* mod, int* fftlist) {
    int i = 0;
    while (mod[i] != '\0') {
        if (!isupper(mod[i])) {
            printf("Error while parsing modification formula:%s\n", mod);
            fflush(stdout);
            return;
        }

        char symbol[3] = {0};
        symbol[0] = mod[i++];
        if (islower(mod[i])) {
            symbol[1] = mod[i++];
        }

        int count = 0;
        while (isdigit(mod[i])) {
            count = count * 10 + (mod[i++] - '0');
        }

        if (count == 0) count = 1;

        for (int j = 0; j < num_simp_elements; j++) {
            if (strcmp(symbol, pep_encoding_elements[j]) == 0) {
                fftlist[j] += count;
                break;
            }
        }
    }
}

//fft
int pep_seq_to_fftlist(const char* sequence, int* fftlist)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 0; // Carbon
    fftlist[1] = 2; // Hydrogen
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 1; // Oxygen
    fftlist[4] = 0; // Sulfur

    int length = strlen(sequence);
    int len = 0;

    int in_mod = 0;
    char curr_mod[100] = {0};
    int mod_index = 0;

    for (int i = 0; i < length; i++) {
        if (in_mod == 0) {
            if (sequence[i] == '[') {
                in_mod = 1;
                continue;
            }

            //Handle the case where the character corresponds to an amino acid.
            int aaindex = aa_to_index(sequence[i]);
            if (aaindex != -1) {
                for (int j = 0;j<num_simp_elements;j++) {
                    fftlist[j] += aa_vectors[aaindex][j];
                }
            }
            len++;
        }
        else {
            if (sequence[i] == ']') {
                in_mod = 0;
                curr_mod[mod_index] = '\0';

                add_mod_to_fftlist(curr_mod, fftlist);

                mod_index = 0;
                memset(curr_mod, 0, sizeof(curr_mod));
            }
            else {
                curr_mod[mod_index] = sequence[i];
                mod_index++;
            }
        }
    }
    return len;
}


int nn_pep_mass_to_isolen(const float mass) {
    if (mass < 1200) {
        return 8;
    }

    if (mass < 11000) {
        return 32;
    }
    if (mass < 55000) {
        return 64;
    }

    if (mass < 120000) {
        return 128;
    }
    return -1;
}


int fft_pep_mass_to_isolen(const float mass) {
    if (mass < 50000) {
        return 64;
    }
    if (mass < 120000) {
        return 128;
    }
    return 1024;
}

//fft
void pep_mass_to_fftlist(const float mass, int* fftlist)
{
    int num = (int)round(mass/massavgine);
    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(num * avgine[i]);
    }
}


int get_pep_isolen_from_seq(const char* seq) {
    int aas = pep_seq_to_aacount(seq);

    if (aas < 50){return 16;}
    if (aas < 300){return 64;}
    return 128;
}


//fft
float fft_pep_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    int* formulalist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return -1.0f;
    }
    pep_seq_to_fftlist(sequence, formulalist);

    int fft_isolen;

    int aas = pep_seq_to_aacount(sequence);

    if (aas > 1000) {
        printf("Sequence contains too many amino acids (>1000).");
        free(formulalist);
        return -1.0f;
    }


    if (aas <= 300) {
        fft_isolen = 64;
    }
    else if (aas <= 1000) {
        fft_isolen = 128;
    }
    else {
        fft_isolen = 512;
    }


    float* fft_isodist = (float*)calloc(fft_isolen, sizeof(float));
    if (fft_isodist == NULL) {
        free(formulalist);
        return -1.0f;
    }

    float maxval = fft_list_to_dist(formulalist, fft_isolen, fft_isodist);
    if (maxval < 0.0f) {
        free(fft_isodist);
        free(formulalist);
        return -1.0f;
    }

    if (fft_isolen < isolen) {
        for (int i = fft_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(fft_isodist);
    free(formulalist);

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }

    return maxval;
}


//fft
float fft_pep_mass_to_dist(const float mass, float *isodist, const int isolen, const int offset)
{
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    int* fftlist = (int*)calloc(5, sizeof(int));
    if (fftlist == NULL) {
        return -1.0f;
    }
    pep_mass_to_fftlist(mass, fftlist);

    int fft_isolen = fft_pep_mass_to_isolen(mass);

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

    if (fft_isolen < isolen) {
        for (int i = fft_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(fft_isodist);
    free(fftlist);

    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }

    return max_val;
}


//nn
float nn_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        return -1.0f;
    }

    mass_to_vector(mass, vector);

    int nn_isolen = nn_pep_mass_to_isolen(mass);

    if (nn_isolen == -1) {
        printf("Error: Mass outside of allowed NN mass range: %f\n", mass);
        free(vector);
        return -1.0f;
    }

    float* nn_isodist = (float*)calloc(nn_isolen, sizeof(float));

    struct IsoGenWeights weights = SetupWeights(5, nn_isolen);
    if (nn_isolen == 8){ weights = LoadWeights(weights, isogenmass_model_8_bin); }
    else if ( nn_isolen == 32 ){ weights = LoadWeights(weights, isogenmass_model_32_bin); }
    else if ( nn_isolen == 64 ){ weights = LoadWeights(weights, isogenmass_model_64_bin); }
    else { weights = LoadWeights(weights, isogenmass_model_128_bin); }

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

float nn_pep_mass_to_dist_custom(const float mass, float* isodist, const int isolen, const int offset,
                                 const char* model_path) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0) {
        return -1.0f;
    }

    float vector[5] = {0};
    mass_to_vector(mass, vector);
    return isogen_model_to_dist_from_file(vector, 5, isodist, isolen, offset, model_path);
}



float brain_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    if (prepare_isodist_output(isodist, isolen, offset) != 0 || !isfinite(mass) || mass <= 0.0f) {
        return -1.0f;
    }
    int* brain_list = (int*)calloc(5, sizeof(int));
    if (brain_list == NULL) {
        return -1.0f;
    }
    pep_mass_to_fftlist(mass, brain_list);

    int brain_isolen = fft_pep_mass_to_isolen(mass);

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


float brain_pep_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset) {
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
    int len = pep_seq_to_fftlist(sequence, formulalist);
    if (len < 1) {
        free(formulalist);
        return -1.0f;
    }

    int brain_isolen;

    if (len <= 300) {
        brain_isolen = 64;
    }
    else if (len <= 1000) {
        brain_isolen = 128;
    }
    else {
        brain_isolen = 1024;
    }

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


#define MAX_LINE_LENGTH 1024
#define INITIAL_CAPACITY 10

typedef struct {
    double mass;
    char* sequence;  // dynamically allocated string
} ProteinEntry;

typedef struct {
    ProteinEntry* entries;
    size_t size;
    size_t capacity;
} ProteinList;

static int parse_protein_entry_line(
    const char* line,
    double* mass,
    char* sequence,
    const size_t sequence_capacity)
{
    char* mass_end = NULL;
    const char* cursor = line;

    while (*cursor != '\0' && isspace((unsigned char)*cursor)) {
        cursor++;
    }
    if (*cursor == '\0' || *cursor == '#') {
        return 0;
    }

    *mass = strtod(cursor, &mass_end);
    if (mass_end == cursor) {
        return 0;
    }

    cursor = mass_end;
    while (*cursor != '\0' && isspace((unsigned char)*cursor)) {
        cursor++;
    }
    if (*cursor == '\0') {
        return 0;
    }

    const size_t sequence_length = strcspn(cursor, " \t\r\n");
    if (sequence_length == 0 || sequence_length >= sequence_capacity) {
        fprintf(stderr, "Skipping malformed sequence entry\n");
        return 0;
    }

    memcpy(sequence, cursor, sequence_length);
    sequence[sequence_length] = '\0';
    return 1;
}

void free_protein_list(ProteinList* list) {
    for (size_t i = 0; i < list->size; ++i) {
        free(list->entries[i].sequence);
    }
    free(list->entries);
    list->entries = NULL;
    list->size = 0;
    list->capacity = 0;
}

ProteinList read_masses_seqs_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        ProteinList empty = {NULL, 0, 0};
        return empty;
    }

    ProteinList list;
    list.size = 0;
    list.capacity = INITIAL_CAPACITY;
    list.entries = malloc(list.capacity * sizeof(ProteinEntry));
    if (!list.entries) {
        fclose(file);
        fprintf(stderr, "Memory allocation failed!\n");
        ProteinList empty = {NULL, 0, 0};
        return empty;
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file)) {
        double mass;
        char seq_buf[MAX_LINE_LENGTH];
        const int parsed = parse_protein_entry_line(
            line,
            &mass,
            seq_buf,
            sizeof(seq_buf)
        );

        if (parsed <= 0) {
            continue;  // skip malformed lines
        }

        if (list.size >= list.capacity) {
            list.capacity *= 2;
            ProteinEntry* temp = realloc(list.entries, list.capacity * sizeof(ProteinEntry));
            if (!temp) {
                fprintf(stderr, "Reallocation failed!\n");
                break;
            }
            list.entries = temp;
        }

        // Allocate and copy the sequence
        list.entries[list.size].mass = mass;
        const size_t sequence_length = strlen(seq_buf) + 1;
        list.entries[list.size].sequence = malloc(sequence_length);
        if (!list.entries[list.size].sequence) {
            fprintf(stderr, "String allocation failed!\n");
            break;
        }
        memcpy(list.entries[list.size].sequence, seq_buf, sequence_length);
        list.size++;
    }

    fclose(file);
    return list;
}



extern void run_file(const char* filename) {
    ProteinList list = read_masses_seqs_file(filename);

    int isolen = 64;

    for (size_t i = 0; i < list.size; i++) {
        float* isodist = (float*)calloc(isolen, sizeof(float));
        if (isodist == NULL) {
            fprintf(stderr, "Error allocating memory for isodist\n");
            break;
        }
        nn_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        fft_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        nn_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        fft_pep_seq_to_dist(list.entries[i].sequence, isodist, isolen, 0);
        nn_pep_seq_to_dist(list.entries[i].sequence, isodist, isolen, 0);
        free(isodist);
    }
    free_protein_list(&list);
    return;
}
