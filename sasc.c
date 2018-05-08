#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mt19937ar.h"
#include <time.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include "tree.h"
#include "utils.h"
#include "sastep.h"
#include "vector.h"

int MAX_ID_TREE = 0;

int main (int argc, char **argv)
{
    // Get CLI arguments
    args_t *arguments = get_arguments(argc, argv);

    int N = arguments->n;
    int M = arguments->m;
    int K = arguments->k;

    double ALPHA = arguments->alpha;
    double BETA = arguments->beta;

    int MAX_LOSSES = arguments->max_del;

    char OUT_PATH[255];
    sprintf(OUT_PATH, "%s_mlt.gv", remove_extension(arguments->infile));

    printf("Starting SASC.\n");

    // Set random seed
    srand((unsigned)time(NULL));
    
    // Create mutation names
    char MUT_NAMES[M][50];
    if (strcmp(arguments->mut_file, "NULL") == 0) {
        for (int i = 0; i < M; i++) {
            sprintf(MUT_NAMES[i], "%d", i+1);
        }
    } else {
        FILE *fp;
        char buff[50];

        fp = fopen(arguments->mut_file, "r");
        int i = 0;
        if (fp == NULL){
            fprintf(stderr, "ERROR: Cannot open MUTFILE.\n");
            exit(EXIT_FAILURE);
        }
            
        while (fgets(buff, sizeof buff, fp) != NULL) {
            // printf("-%s-\n", buff);
            buff[strcspn(buff, "\n")] = '\0';
            strcpy(MUT_NAMES[i], buff);
            i++;
        }
        fclose(fp);

        if (i != M) {
            fprintf(stderr, "ERROR: Dimesion of mutations names and MUTATIONS are different.\n");
            exit(EXIT_FAILURE);
        }
            
    }

    //Itialize MT19937

    unsigned long init[10], length=10;
    FILE *f;
    f = fopen("/dev/urandom", "r");

    for (int i = 0; i < length; i++){
        unsigned long randval;
        fread(&randval, sizeof(randval), 1, f);
        init[i] = randval;
    }
    fclose(f);

    init_by_array(init, length);

    // Import INFILE
    int **INPUT_MATRIX;
    INPUT_MATRIX = malloc(N * sizeof *INPUT_MATRIX);
    for (int i = 0; i < N; i++) {
        INPUT_MATRIX[i] = malloc(M * sizeof *INPUT_MATRIX[i]);
    }

    import_input(INPUT_MATRIX, N, M, arguments->infile);


    double START_TEMP = 100000.0;
    double COOLING_RATE = 0.00001;    
    double MIN_TEMP = 0.0001;
    // These three are for debug purposes    
    START_TEMP = 10000.0;
    COOLING_RATE = 0.01;
    MIN_TEMP = 0.001;

    int REPETITIONS = arguments->repetitions;
    node_t *best_tree = NULL;
    double best_loglike = -DBL_MAX;
    int best_sigma[N]; for (int i = 0; i < N; i++) { best_sigma[i] = 0; }
    vector best_tree_vec;
    vector_init(&best_tree_vec);
    vector best_losses_vec;
    vector_init(&best_losses_vec);

    for (int r = 0; r < REPETITIONS; r++) {
        printf("Iteration: %d\n", r+1);

        // Generate a RANDOM BTREE
        vector ml_tree_vec;
        vector_init(&ml_tree_vec);

        node_t *root = node_new("germline", -1, 0);
        vector_add(&ml_tree_vec, root);

        int rantree[M];
        for (int i = 0; i < M; i++) {
            rantree[i] = i;
        }
        shuffle(rantree, M);

        int app_node = 0;
        for (int i = 0; i < M; i++) {
            node_t *cnode1 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
            vector_add(&ml_tree_vec, cnode1);
            node_append(vector_get(&ml_tree_vec, app_node), cnode1);
            i += 1;

            if (i < M){
                node_t *cnode2 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&ml_tree_vec));
                vector_add(&ml_tree_vec, cnode2);
                node_append(vector_get(&ml_tree_vec, app_node), cnode2);
            }
            app_node += 1;
        }

        // Generate SIGMA
        int SIGMA[N];
        for (int i = 0; i < N; i++) {
            SIGMA[i] = 0;
        }


        // get log-likelihood
        double lh = greedy_tree_loglikelihood(root, ml_tree_vec, SIGMA, INPUT_MATRIX, N, M, ALPHA, BETA);
        // for (int i = 0; i < N; i++) { printf("%d ", SIGMA[i]); } printf("\n");
        // printf("Start log-like: %lf\n", lh);
        node_t *ml_tree = anneal(root, ml_tree_vec, N, M, K, ALPHA, BETA, INPUT_MATRIX, START_TEMP, COOLING_RATE, MIN_TEMP, MAX_LOSSES);
        
        vector_free(&ml_tree_vec);
        vector_init(&ml_tree_vec);
        vector ml_losses_vec;
        vector_init(&ml_losses_vec);
        
        ml_tree = treecpy(ml_tree, &ml_tree_vec, &ml_losses_vec, N);
        
        double current_lh = greedy_tree_loglikelihood(ml_tree, ml_tree_vec, SIGMA, INPUT_MATRIX, N, M, ALPHA, BETA);
        // printf("Maximum log-likelihood: %lf\n", current_lh);

        if (current_lh > best_loglike) {
            if (best_tree != NULL)
                destroy_tree(best_tree);

            vector_free(&best_tree_vec);
            vector_free(&best_losses_vec);

            vector_init(&best_tree_vec);
            vector_init(&best_losses_vec);

            best_tree = treecpy(ml_tree, &best_tree_vec, &best_losses_vec, N);
        }

        vector_free(&ml_tree_vec);
        destroy_tree(root);
        destroy_tree(ml_tree);
        
    }
   
    double best_calculated_likelihood = greedy_tree_loglikelihood(best_tree, best_tree_vec, best_sigma, INPUT_MATRIX, N, M, ALPHA, BETA);    
    // printf("Maximum likelihood found: %f\n", best_calculated_likelihood);

    if (arguments->print_leaves == 1) {
        fprint_tree_leaves(best_tree, &best_tree_vec, best_sigma, N, OUT_PATH);
    } else {
        fprint_tree(best_tree, OUT_PATH);
    }
    
    printf("\nMaximum likelihood Tree found: %f\n", best_calculated_likelihood);
    print_tree(best_tree);    
    
    printf("Best cell assigment:\n");
    for (int i=0; i<N;i++) {
        printf("%d ", best_sigma[i]);
    }
    printf("\n");

    if (arguments->print_expected == 1) {
        char OUT_MATRIX[255];
        sprintf(OUT_MATRIX, "%s_out.txt", remove_extension(arguments->infile));

        FILE *fpo;
        fpo = fopen(OUT_MATRIX, "w+");
        int gtpo[M];

        for (int i=0; i<N;i++) {
            for (int j = 0; j < M; j++) { gtpo[j] = 0; }

            node_t *node = vector_get(&best_tree_vec, best_sigma[i]);
            if (node == NULL) printf("i: %d --- sigma: %d", i, best_sigma[i]);
            assert(node != NULL);
            get_genotype_profile(vector_get(&best_tree_vec, best_sigma[i]), gtpo);

            for (int j = 0; j < M; j++) {
                fprintf(fpo, "%d ", gtpo[j]);
            }
            fprintf(fpo, "\n");

        }
        fclose(fpo);
    }


    return 0;
}