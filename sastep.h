#include "vector.h"

int random_assignment(int MAX);

// int prune_regraft(node_t *prune, node_t *regraft);
// int add_back_mutation(node_t *node, node_t *tree[], int max_tree_id, int m);

double tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double alpha, double beta);

node_t * anneal(node_t *root, int sigma[], vector tree_vec, int n, int m, int k, double alpha, double beta,  int **inmatrix, double start_temp, double cooling_rate, double min_temp, double *best_lh, int MAX_LOSSES);