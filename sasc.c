#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mt19937ar.h"
#include <time.h>
#include <assert.h>
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

    char OUT_PATH[255];
    sprintf(OUT_PATH, "%s_mlt.gv", remove_extension(arguments->infile));

  /* Print argument values */
    printf("n: %d\n", N);
    printf("m: %d\n", M);
    printf("k: %d\n", K);
    printf("a: %f\n", ALPHA);
    printf("b: %f\n", BETA);
    printf("i: %s\n", arguments->infile);
    printf("e: %s\n", arguments->mut_file);

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

    // Generate a RANDOM BTREE
    node_t *TREE[K*M]; // delete this
    vector tree_nodes;
    vector_init(&tree_nodes);

    node_t *root = node_new("germline", -1, 0);
    vector_add(&tree_nodes, root);

    int rantree[M];
    for (int i = 0; i < M; i++) {
        rantree[i] = i;
    }
    shuffle(rantree, M);

    int app_node = 0;
    for (int i = 0; i < M; i++) {
        node_t *cnode1 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&tree_nodes));
        vector_add(&tree_nodes, cnode1);
        node_append(vector_get(&tree_nodes, app_node), cnode1);
        i += 1;

        if (i < M){
            node_t *cnode2 = node_new(MUT_NAMES[rantree[i]], rantree[i], vector_total(&tree_nodes));
            vector_add(&tree_nodes, cnode2);
            node_append(vector_get(&tree_nodes, app_node), cnode2);
        }
        app_node += 1;
    }

//    printf("m: %d, v_tot: %d\n", M, vector_total(&tree_nodes));

    for (int i = 0; i < vector_total(&tree_nodes); i++) {

        node_t *n = vector_get(&tree_nodes, i);
        assert(n->id == i);
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

    // Generate SIGMA
    int SIGMA[N];



    for (int i = 0; i < N; i++) {
        SIGMA[i] = random_assignment(M);
    }
    
    // fprint_tree(root, OUT_PATH);
    // fprint_tree_leaves(root, TREE, SIGMA, N, OUT_PATH);

//     for (int i=0; i<N;i++) {
//         printf("%d,", SIGMA[i]);
//     }
//     printf("\n");

    // Import INFILE
    int **INPUT_MATRIX;
    INPUT_MATRIX = malloc(N * sizeof *INPUT_MATRIX);
    for (int i = 0; i < N; i++) {
        INPUT_MATRIX[i] = malloc(M * sizeof *INPUT_MATRIX[i]);
    }

    import_input(INPUT_MATRIX, N, M, arguments->infile);

//     print_tree(root);

    // get log-likelihood
    double lh = tree_loglikelihood(root, tree_nodes, SIGMA, INPUT_MATRIX, N, M, ALPHA, BETA);
    printf("Start log-like: %lf\n", lh);

    double START_TEMP = 100000.0;
    double COOLING_RATE = 0.00001;
    double MIN_TEMP = 0.000001;
    // int REPETITIONS = 1;

    // node_t *ml_tree = root;

    // for (int r = 0; r < REPETITIONS; r++) {
    //     printf("Iteration: %d\n", r+1);
    node_t *ml_tree = anneal(root, SIGMA, tree_nodes, N, M, K, ALPHA, BETA, INPUT_MATRIX, START_TEMP, COOLING_RATE, MIN_TEMP);
    // }
   
    
//    fprint_tree_leaves(ml_tree, TREE, SIGMA, N, OUT_PATH);
    print_tree(ml_tree);
    
    printf("Cell assigment:\n");
    for (int i=0; i<N;i++) {
        printf("%d ", SIGMA[i]);
    }
    printf("\n");

//    int gtp[M];
//    for (int j = 0; j < M; j++) {
//            gtp[j] = 0;
//        }
//
//        get_genotype_profile(vector_get(&tree_nodes, SIGMA[32]), gtp);
//
//    printf("Genotype profile of cell32:\n");
//    for (int j = 0; j < M; j++) {
//        printf("%d ", gtp[j]);
//    }
//    printf("\n");

    return 0;
}