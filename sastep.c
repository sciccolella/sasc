/*
    MIT License

    Copyright (c) 2017-2019 Simone Ciccolella

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "vector.h"
#include "tree.h"
#include "sastep.h"
#include "mt19937ar.h"
#include <omp.h>

#ifdef NDEBUG
#include <assert.h>
#else
#define assert(ignore)((void) 0)
#endif

int
random_assignment(int MAX){
    unsigned long randval = genrand_int32();
    return (int) (randval % (MAX + 1));
}

void
check_subtree_losses(node_t *node, vector *tree_vec, vector *loss_vec, int *k_loss, int *sigma, int n) {
    if (node == NULL) return;

    check_subtree_losses(node->first_child, tree_vec, loss_vec, k_loss, sigma, n);
    check_subtree_losses(node->next_sibling, tree_vec, loss_vec, k_loss, sigma, n);

    if (node->loss == 1){
        bool valid = is_loss_valid(node);
        bool lost = is_already_lost(node, node->mut_index);

        if (valid == false || lost == true) {
            node_delete(node, tree_vec, loss_vec, k_loss, sigma, n);
        }

    }
}

int
prune_regraft(node_t *prune, node_t *regraft, node_t *root, int monoclonal) {
    if (is_ancestor(regraft, prune) == true)
        return 1;
    if (regraft->parent == prune)
        return 1;
    if (prune->parent == NULL)
        return 1;
    if (prune == regraft)
        return 1;

    if (monoclonal == 1 && regraft->parent == NULL)
        return 1;

    node_detach(prune);
    node_append(regraft, prune);
    
    return 0;
}

int
add_back_mutation(node_t *node, vector *tree_vec, int m, int k, int *k_loss, vector *losses_vec, int MAX_LOSSES) {
    node_t *par = node->parent;
    if (par == NULL || par->parent == NULL)
        return 1;

    if (vector_total(losses_vec) >= MAX_LOSSES)
        return 1;

    // Walk to root to select possible candidates for deletion
    node_t *candidates[m*k+1];
    for (int i =0; i<k*m+1; i++) { candidates[i] = NULL; }

    int tot = 0;
    while (par != NULL) {
        if (par->loss == 0) {
            candidates[tot] = par;
            tot++;
        }
        par = par->parent;
    }

    int rand_del = random_assignment(tot-2) + 1;
    if (candidates[rand_del] == NULL)
        return 1;
    if (candidates[rand_del]->mut_index == -1)
        return 1;
    
    char label[255];
    strcpy(label, candidates[rand_del]->label);
    if (k_loss[candidates[rand_del]->mut_index] >= k)
        return 1;
    if (is_already_lost(node, candidates[rand_del]->mut_index) == true)
        return 1;


    node_t *node_del = node_new(label, candidates[rand_del]->mut_index, vector_total(tree_vec));
    node_del->loss = 1;

    vector_add(tree_vec, node_del);
    vector_add(losses_vec, node_del);

    k_loss[node_del->mut_index] += 1;

    par = node->parent;
    node_detach(node);
    node_append(par, node_del);
    node_append(node_del, node);

    return 0;
}

double
greedy_tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double* alpha,
        double beta, double *gammas, int *k_loss, int CORES) {

    int node_max = vector_total(&tree_vec);

    int* nodes_genotypes = calloc(node_max * m, sizeof(int));

    for (int i = 0; i < node_max; i++) {

        node_t *node = vector_get(&tree_vec, i);

        if (node == NULL) {
            for (int j = 0; j < m; j++) {
                nodes_genotypes[i * m + j] = 3;
            }
        } else {
            get_genotype_profile(node, nodes_genotypes + i * m);
        }
    }

    double loss_weight = 0;
    for (int j = 0; j < m; j++) {
        loss_weight += k_loss[j] * log(gammas[j]);
    }

    double maximum_likelihood = 0 - loss_weight;

    double like_00 = log(1-beta);
    double like_10 = log(beta);
    double *like_01 = malloc(m * sizeof(double));
    double *like_11 = malloc(m * sizeof(double));
    double like_2 = 0;

    for (int j = 0; j < m; j++) {
        like_01[j] = log(alpha[j]);
        like_11[j] = log(1 - alpha[j]);
    }

    #pragma omp parallel for reduction(+:maximum_likelihood) num_threads(CORES)
    for (int i = 0; i < n; i++) {
        int best_sigma = -1;
        double best_lh = -DBL_MAX;

        for (int node = 0; node < node_max; node++) {
        // for (node = 0; node < node_max; node++) {
            if (nodes_genotypes[node  * m + 0] != 3) {
                double lh = 0;

                for (int j = 0; j < m; j++) {
                    double p = 1;
                    // double p = prob(inmatrix[i][j], nodes_genotypes[node * m + j], alpha[j], beta);
                    // int I = inmatrix[i][j];
                    // int E = nodes_genotypes[node * m + j];
                    if (inmatrix[i][j] == 0 && nodes_genotypes[node * m + j] == 0) {
                        p = like_00;
                    } else if (inmatrix[i][j] == 0 && nodes_genotypes[node * m + j] == 1) {
                        p = like_01[j];
                    } else if (inmatrix[i][j] == 1 && nodes_genotypes[node * m + j] == 0) {
                        p = like_10;
                    } else if (inmatrix[i][j] == 1 && nodes_genotypes[node * m + j] == 1) {
                        p = like_11[j];
                    } else if (inmatrix[i][j] == 2) {
                        p = like_2;
                    }
                    lh += p;
                }

                if (lh > best_lh) {
                    best_sigma = node;
                    best_lh = lh;
                }
            }
        }
        sigma[i] = best_sigma;
        maximum_likelihood += best_lh;
    }

    free(nodes_genotypes);
    free(like_01);
    free(like_11);
    return maximum_likelihood;
}

double
gen_gaussian(double mu, double sigma) {
    static const double two_pi = 2.0*3.14159265358979323846;
    double u1 = genrand_real3();
    double u2 = genrand_real3();
    double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    return fabs(z0 * sigma + mu);

}

void el_commit(elpar_t* el_params, double* beta) {
    for (int i = 0; i < el_params->M; i++) {
        el_params->ALPHAS[i] = el_params->a_xs[i];
        el_params->GAMMAS[i] = el_params->g_xs[i];
    }
    *beta = el_params->b_x;
    el_params->changed = 0;

}

void el_discard(elpar_t* el_params, double beta) {
    for (int i = 0; i < el_params->M; i++) {
        el_params->a_xs[i] = el_params->ALPHAS[i];
        el_params->g_xs[i] = el_params->GAMMAS[i];
    }
    el_params->b_x = beta;
    el_params->changed = 0;
}

void
learn_a(elpar_t* el_params) {
    double x = 0;
    if (el_params->single_alpha == 1) {
        x = gen_gaussian(el_params->a_mu[0], el_params->a_variance);
        for (int i = 0; i < el_params->M; i++) {
            el_params->a_xs[i] = x;
        }
    } else {
        for (int i = 0; i < el_params->M; i++) {
            x = gen_gaussian(el_params->a_mu[i], el_params->a_variance);
            el_params->a_xs[i] = x;
        }
    }
}

void
learn_b(elpar_t* el_params) {
    double x = 0;
    x = gen_gaussian(el_params->b_mu, el_params->b_variance);
    el_params->b_x = x;
}

void
learn_g(elpar_t* el_params) {
    double x = 0;
    if (el_params->single_gamma == 1) {
        x = gen_gaussian(el_params->g_mu[0], el_params->g_variance);
        for (int i = 0; i < el_params->M; i++) {
            el_params->g_xs[i] = x;
        }
    } else {
        for (int i = 0; i < el_params->M; i++) {
            x = gen_gaussian(el_params->g_mu[i], el_params->g_variance);
            el_params->g_xs[i] = x;
        }
    }
}

void
params_learning(elpar_t* el_params) {
    double choice = genrand_real1();

    if (el_params->a_variance > 0 && el_params->b_variance > 0 && el_params->g_variance > 0) {
        if (choice < 0.33) {
            learn_a(el_params);
        } else if (choice < 0.66) {
            learn_b(el_params);
        } else {
            learn_g(el_params);
        }
    } else if (el_params->a_variance > 0 && el_params->b_variance > 0) {
        if (choice < 0.5) {
            learn_a(el_params);
        } else {
            learn_b(el_params);
        }
    } else if (el_params->a_variance > 0 && el_params->g_variance > 0) {
        if (choice < 0.5) {
            learn_a(el_params);
        } else {
            learn_g(el_params);
        }
    } else if (el_params->b_variance > 0 && el_params->g_variance > 0) {
        if (choice < 0.5) {
            learn_b(el_params);
        } else {
            learn_g(el_params);
        }
    } else if (el_params->a_variance > 0) {
        learn_a(el_params);
    } else if (el_params->b_variance > 0) {
        learn_b(el_params);
    } else if (el_params->g_variance > 0) {
        learn_g(el_params);
    }

    el_params->changed = 1;
}

void
neighbor(node_t *root, vector *tree_vec, int *sigma, int m, int n, int k, vector *loss_vec, int *k_loss, int MAX_LOSSES,
            elpar_t* el_params, int monoclonal) {

    double el = genrand_real1();
    if (el < 0.1 && (el_params->a_variance > 0 || el_params->b_variance > 0 || el_params->g_variance > 0)) {
        params_learning(el_params);
    } else {
        double move = genrand_real1();
        if (move < 0.25 && k > 0) {
            // Add back-mutation
            int bm_res = 1;
            node_t *node_res = NULL;

            int ip = random_assignment(vector_total(tree_vec) - 1);
            node_res = vector_get(tree_vec, ip);

            bm_res = add_back_mutation(node_res, tree_vec, m, k, k_loss, loss_vec, MAX_LOSSES);

            if (bm_res == 0) {
                check_subtree_losses(node_res, tree_vec, loss_vec, k_loss, sigma, n);
            }
        } else if (move < 0.50 && k > 0) {
            // Delete a mutation
            node_t *node_res = NULL;
            if (vector_total(loss_vec) == 0)
                return;

            int node_max = vector_total(loss_vec) - 1;
            assert(node_max >= 0);
            int ip = random_assignment(node_max);
            node_res = vector_get(loss_vec, ip);

            node_delete(node_res, tree_vec, loss_vec, k_loss, sigma, n);

        } else if ((move < 0.75 && k > 0) || (move < 0.50 && k == 0)) {
            // switch nodes
            node_t *u = NULL;
            while (u == NULL || u->parent == NULL || u->loss == 1) {
                int node_max = vector_total(tree_vec) - 1;
                assert(node_max > 0);
                int ip = random_assignment(node_max);
                u = vector_get(tree_vec, ip);

            }

            node_t *v = NULL;
            while (v == NULL || v->parent == NULL || v->loss == 1 || v->id == u->id) {
                int node_max = vector_total(tree_vec) - 1;
                assert(node_max > 0);
                int ig = random_assignment(node_max);
                v = vector_get(tree_vec, ig);
            }

            int mut_tmp;
            char label_tmp[255];

            mut_tmp = u->mut_index;
            strcpy(label_tmp, u->label);

            u->mut_index = v->mut_index;
            strcpy(u->label, v->label);

            v->mut_index = mut_tmp;
            strcpy(v->label, label_tmp);

            check_subtree_losses(u, tree_vec, loss_vec, k_loss, sigma, n);
            check_subtree_losses(v, tree_vec, loss_vec, k_loss, sigma, n);
        } else {
            // Prune-regraft two random nodes
            int pr_res = 1;
            node_t *prune_res = NULL;
            while (pr_res != 0) {
                node_t *prune = NULL;
                while (prune == NULL || prune->parent == NULL) {
                    int node_max = vector_total(tree_vec) - 1;
                    assert(node_max > 0);
                    int ip = random_assignment(node_max);
                    prune = vector_get(tree_vec, ip);

                }

                node_t *graft = NULL;
                while (graft == NULL) {
                    int node_max = vector_total(tree_vec) - 1;
                    assert(node_max > 0);
                    int ig = random_assignment(node_max);
                    graft = vector_get(tree_vec, ig);
                }
                pr_res = prune_regraft(prune, graft, root, monoclonal);
                prune_res = prune;
            }
            check_subtree_losses(prune_res, tree_vec, loss_vec, k_loss, sigma, n);
        }
    }
}

elpar_t* set_el_params(int single_a, int m, double* ALPHAS, double* a_mu, double a_variance, double* a_xs,
                       double* beta, double b_mu, double b_variance, double* GAMMAS, double* g_mu, double g_variance, double* g_xs, int single_g){

    elpar_t *params = malloc(sizeof(elpar_t));
    params->single_alpha = single_a;
    params->M = m;
    params->changed = 0;

    params->ALPHAS = ALPHAS;
    params->a_mu = a_mu;
    params->a_variance = a_variance;
    params->a_xs = a_xs;

    params->BETA = beta;
    params->b_mu = b_mu;
    params->b_variance = b_variance;
    params->b_x = b_mu;

    params->single_gamma = single_g;

    params->GAMMAS = GAMMAS;
    params->g_mu = g_mu;
    params->g_variance = g_variance;
    params->g_xs = g_xs;

    return params;
}

double
accept_prob(double old_lh, double new_lh, double current_temp) {
    if (new_lh > old_lh) {
        return 1.0;
    } else {
        double a = exp((new_lh - old_lh) / current_temp);
        return a;
    }
}

node_t *
anneal(node_t *root, vector tree_vec, int n, int m, int k, double* alpha, double beta,  int **inmatrix,
       double start_temp, double cooling_rate, double min_temp, int MAX_LOSSES, elpar_t* el_params, double* gamma, int *Cj, int MONOCLONAL, int CORES) {

    double current_temp = start_temp;
    double current_cooling_rate = cooling_rate;

    // Vector of node indices
    vector current_tree_vec;
    vector_init(&current_tree_vec);

    // Vector of back-mutation node indices
    vector current_losses_vec;
    vector_init(&current_losses_vec);

    // Current count of losses per mutation to ensure them to be < k
    int current_kloss[m]; for (int i =0; i< m; i ++) { current_kloss[i] = 0; }

    // Current assignment
    int current_sigma[n];
    for (int i = 0; i < n; i++) { current_sigma[i] = 0; }

    node_t *current_root = treecpy(root, &current_tree_vec, &current_losses_vec, n);

    unsigned int step = 0;

    double current_lh = greedy_tree_loglikelihood(current_root, tree_vec, current_sigma, inmatrix, n, m, alpha, beta, gamma, current_kloss, CORES);

    printf("Step\t\t\tLog-likelihood\t\t\tTemperature\n");
    while (current_temp > min_temp) {
        // Create a modifiable copy
        vector copy_tree_vec;
        vector_init(&copy_tree_vec);

        vector copy_losses_vec;
        vector_init(&copy_losses_vec);
        int copy_kloss[m];
        for (int i = 0; i < m; i++) { copy_kloss[i] = current_kloss[i]; }

        int copy_sigma[n];
        for (int i = 0; i < n; i++) { copy_sigma[i] = current_sigma[i]; }

        node_t *copy_root = treecpy(current_root, &copy_tree_vec, &copy_losses_vec, n);


        for (int i = 0; i < vector_total(&copy_tree_vec); i++) {
            node_t *n = vector_get(&copy_tree_vec, i);
            assert(n->id == i);
        }
        assert(vector_total(&copy_tree_vec) == vector_total(&current_tree_vec));
        assert(vector_total(&copy_losses_vec) == vector_total(&current_losses_vec));

        neighbor(copy_root, &copy_tree_vec, copy_sigma, m, n, k, &copy_losses_vec, copy_kloss, MAX_LOSSES, el_params, MONOCLONAL);

        double new_lh = 0;
        if (el_params->changed == 1) {
            // TODO: change gamma to g_x
            new_lh = greedy_tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma, inmatrix, n, m, el_params->a_xs, el_params->b_x, el_params->g_xs, copy_kloss, CORES);
        } else {
            new_lh = greedy_tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma, inmatrix, n, m, alpha, beta, gamma, copy_kloss, CORES);
        }


        double acceptance = accept_prob(current_lh, new_lh, current_temp);
        double p = genrand_real1();

        if (acceptance > p) {

            if (el_params->changed == 1) {
                el_commit(el_params, &beta);
            }
            // Destroy current solution
            destroy_tree(current_root);

            vector_free(&current_tree_vec);
            vector_init(&current_tree_vec);

            vector_free(&current_losses_vec);
            vector_init(&current_losses_vec);

            // Copy new solution to current
            double test_lh = greedy_tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma, inmatrix, n, m, alpha, beta, gamma, copy_kloss, CORES);
            assert(test_lh == new_lh);

            current_lh = new_lh;

            for (int i = 0; i < m; i++) { current_kloss[i] = copy_kloss[i]; }

            current_root = treecpy(copy_root, &current_tree_vec, &current_losses_vec, n);
            double test_clh = greedy_tree_loglikelihood(current_root, current_tree_vec, current_sigma, inmatrix, n, m, alpha, beta, gamma, current_kloss, CORES);
            assert(test_clh == test_lh);


        } else {
            if (el_params->changed == 1) {
                el_discard(el_params, beta);
            }
        }

        // Destroy neighbor
        destroy_tree(copy_root);
        vector_free(&copy_tree_vec);
        vector_free(&copy_losses_vec);


        
        if (current_temp <= 0.25 * start_temp) {
            current_temp *= (1 - (current_cooling_rate * 0.1));
        } else {
            current_temp *= (1 - current_cooling_rate);
        }

        ++step;
        if (step % 1000 == 0 || step == 1) {
            printf("%d\t\t\t%lf\t\t\t%lf\n", step, current_lh, current_temp);
        }

    }


    // Change SIGMA to best solution
    int bs[n];
    for (int i = 0; i < n; i++) { bs[i] = 0; }

    vector_init(&tree_vec);

    vector best_losses_vec;
    vector_init(&best_losses_vec);
    node_t *best_root = treecpy(current_root, &tree_vec, &best_losses_vec, n);

    for (int i = 0; i < m; i++) { Cj[i] = current_kloss[i]; }

    // Return best TREE
    return best_root;
}
