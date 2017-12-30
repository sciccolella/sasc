typedef struct Arguments {
    char infile[255];
    int n;
    int m;
    int k;
    double alpha;
    double beta;
    char mut_file[255];
} args_t;

void shuffle(int *array, int n);
// unsigned int get_randint();
char *remove_extension (char* path);
void import_input(int **input_matrix, int rows, int columns, char *path);
args_t* get_arguments(int cargc, char **cargsv);