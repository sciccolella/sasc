#include "vector.h"

typedef struct el_params {
    int single_alpha;
    int single_gamma;
    int M;

    double* ALPHAS;
    double* a_mu;
    double a_variance;
    double* a_xs;

    double* BETA;
    double b_mu;
    double b_variance;
    double b_x;

    double* GAMMAS;
    double* g_mu;
    double g_variance;
    double* g_xs;

    int changed;
} elpar_t;

int random_assignment(int MAX);

// int prune_regraft(node_t *prune, node_t *regraft);
// int add_back_mutation(node_t *node, node_t *tree[], int max_tree_id, int m);

//double tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double alpha, double beta);
double greedy_tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double* alpha, double beta, double *gammas, int *k_loss);

node_t * anneal(node_t *root, vector tree_vec, int n, int m, int k, double* alpha, double beta,  int **inmatrix,
                double start_temp, double cooling_rate, double min_temp, int MAX_LOSSES, elpar_t* el_params, double* gamma, int *Cj);

elpar_t* set_el_params(int single, int m, double* ALPHAS, double* a_mu, double a_variance, double* a_xs,
                       double* beta, double b_mu, double b_variance, double* GAMMAS, double* g_mu, double g_variance, double* g_xs, int single_g);