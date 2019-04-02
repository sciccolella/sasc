#ifndef UTILS_H
#define UTILS_H

typedef struct Arguments {
    char infile[255];
    int n;
    int m;
    int k;
    int max_del;
    double alpha;
    char alpha_file[255];
    double beta;
    double gamma;
    char gamma_file[255];
    char mut_file[255];
    char cell_file[255];
    int print_leaves;
    int print_expected;
    int monoclonal;
    int repetitions;
    double start_temp;
    double cooling_rate;
    double el_a_variance;
    double el_b_variance;
    double el_g_variance;
} args_t;

void shuffle(int *array, int n);
// unsigned int get_randint();
char *remove_extension (char* path);
void import_input(int **input_matrix, int rows, int columns, char *path);
args_t* get_arguments(int cargc, char **cargsv);

#endif