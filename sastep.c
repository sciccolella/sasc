#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "vector.h"
#include "tree.h"
#include "sastep.h"
#include "mt19937ar.h"

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
//    printf("%s\n", node->label);

    if (node->loss == 1){
        bool valid = is_loss_valid(node);
        bool lost = is_already_lost(node, node->mut_index);

#ifdef DEBUG
        printf("valid: %s\n", valid ? "true" : "false");
        printf("lost: %s\n", lost ? "true" : "false");
#endif

        if (valid == false || lost == true) {
            node_delete(node, tree_vec, loss_vec, k_loss, sigma, n);
        }

    }
}

int
prune_regraft(node_t *prune, node_t *regraft, node_t *root) {
    if (is_ancestor(regraft, prune) == true)
        return 1;
    if (regraft->parent == prune)
        return 1;
    if (prune->parent == NULL)
        return 1;
    if (prune == regraft)
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
    // printf("%d\n", rand_del);
    // printf("%s\n", node->label);
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
prob(int I, int E, double alpha, double beta){
    double p = 0;
    switch(I) {
        case 0:
            switch(E) {
                case 0:
                    p = 1 - beta;
                    break;
                case 1:
                    p = alpha;
                    break;
                default:
                    fprintf(stderr, "ERROR: Unkown value of E (%d).\n", E);
                    exit(EXIT_FAILURE);
            }
            break;
        case 1:
            switch(E) {
                case 0:
                    p = beta;
                    break;
                case 1:
                    p = 1 - alpha;
                    break;
                default:
                    fprintf(stderr, "ERROR: Unkown value of E (%d).\n", E);
                    exit(EXIT_FAILURE);
            }
            break;
        case 2:
            p = 1;
            break;
        default:
            fprintf(stderr, "ERROR: Unkown value of I (%d).\n", I);
            exit(EXIT_FAILURE);
    }

    return p;
}

//double
//tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double alpha, double beta) {
//    double lh = 0;
//    int gtp[m];
//
//    // for (int i=0; i<n;i++) {
//    //     printf("%d,", sigma[i]);
//    // }
//    // printf("\n");
//    // print_tree(root);
//
//    for (int i = 0; i < n; i++) {
//
//
//        for (int j = 0; j < m; j++) {
//            gtp[j] = 0;
//        }
////        printf("cell: %d, sigma: %d\n", i, sigma[i]);
//        node_t *node = vector_get(&tree_vec, sigma[i]);
//        assert(node != NULL);
//        get_genotype_profile(vector_get(&tree_vec, sigma[i]), gtp);
//
////        for (int j = 0; j < m; j++) {
////            printf("%d ", gtp[j]);
////        }
////        printf("\n");
//
//        for (int j = 0; j < m; j++) {
//            double p = prob(inmatrix[i][j], gtp[j], alpha, beta);
//            lh += log(p);
//        }
//    }
//
//    return lh;
//}

double
tree_loglikelihood(node_t *root, vector tree_vec, int *sigma, int **inmatrix, int n, int m, double alpha, double beta) {
    double lh = 0;
    int gtp[m];

    for (int i = 0; i < n; i++) {
        double max_cell_lh = -DBL_MAX;
        int max_cell_node = -1;

        for (int node_id = 0; node_id < vector_total(&tree_vec); node_id++) {

            double curr_lh = 0;

            for (int j = 0; j < m; j++) {
                gtp[j] = 0;
            }

            node_t *node = vector_get(&tree_vec, node_id);
            if (node == NULL) {
                continue;
            }
            get_genotype_profile(node, gtp);

            for (int j = 0; j < m; j++) {
                double p = prob(inmatrix[i][j], gtp[j], alpha, beta);
                curr_lh += log(p);
            }

            if (curr_lh > max_cell_lh) {
                max_cell_lh = curr_lh;
                max_cell_node = node_id;
            }
        }
        sigma[i] = max_cell_node;
        lh += max_cell_lh;
    }

    return lh;
}

void
neighbor(node_t *root, vector *tree_vec, int *sigma, int m, int n, int k, vector *loss_vec, int *k_loss, int MAX_LOSSES, int phase) {
    double move = genrand_real1();
    // printf("%lf,", move);
    if (phase == 2) {
        if (move < 0.15) {
            // Add back-mutation
            int bm_res = 1;
            node_t *node_res = NULL;

            int ip = random_assignment(vector_total(tree_vec) - 1);
            node_res = vector_get(tree_vec, ip);

            bm_res = add_back_mutation(node_res, tree_vec, m, k, k_loss, loss_vec, MAX_LOSSES);

            if (bm_res == 0) {
                check_subtree_losses(node_res, tree_vec, loss_vec, k_loss, sigma, n);
#ifdef DEBUG
                printf("\tADD BACK-MUT: %s (id: %d)\n", node_res->label, node_res->id);
#endif
            }
        } else if (move < 0.30) {
            // Delete a mutation
            node_t *node_res = NULL;
            if (vector_total(loss_vec) == 0)
                return;

            int node_max = vector_total(loss_vec) - 1;
            assert(node_max >= 0);
            int ip = random_assignment(node_max);
            node_res = vector_get(loss_vec, ip);

            node_delete(node_res, tree_vec, loss_vec, k_loss, sigma, n);

#ifdef DEBUG
            printf("\tDELETE BACK-MUT: %s (id: %d)\n", node_res->label, node_res->id);
#endif
//        } else if (move <= 0.65){
//            // Change assigment
//            int rand_cell = random_assignment(n);
//
//            int node_max = vector_total(tree_vec) - 1;
//            assert(node_max > 0);
//            sigma[rand_cell] = random_assignment(node_max);
        } else {
            // switch nodes
            node_t *u = NULL;
            while (u == NULL || u->parent == NULL || u->loss == 1) {
                int node_max = vector_total(tree_vec) - 1;
                assert(node_max > 0);
                int ip = random_assignment(node_max);
                // printf("p:%d\n", ip);
                u = vector_get(tree_vec, ip);

            }

            node_t *v = NULL;
            while (v == NULL || v->parent == NULL || v->loss == 1 || v->id == u->id) {
                int node_max = vector_total(tree_vec) - 1;
                assert(node_max > 0);
                int ig = random_assignment(node_max);
                // printf("g:%d\n", ig);
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
        }
    } else {
//        if (move < 0.5) {
//            // Change assigment
//            int rand_cell = random_assignment(n);
//
//            int node_max = vector_total(tree_vec) - 1;
//            assert(node_max > 0);
//            sigma[rand_cell] = random_assignment(node_max);
//#ifdef DEBUG
//            printf("\tASS: cell: %d, tot: %d, ass: %d\n", rand_cell, node_max, sigma[rand_cell]);
//#endif
//        } else {
            // Prune-regraft two random nodes
            int pr_res = 1;
            node_t *prune_res = NULL;
            while (pr_res != 0) {
                node_t *prune = NULL;
                while (prune == NULL || prune->parent == NULL) {
                    int node_max = vector_total(tree_vec) - 1;
                    assert(node_max > 0);
                    int ip = random_assignment(node_max);
                    // printf("p:%d\n", ip);
                    prune = vector_get(tree_vec, ip);

                }

                node_t *graft = NULL;
                while (graft == NULL) {
                    int node_max = vector_total(tree_vec) - 1;
                    assert(node_max > 0);
                    int ig = random_assignment(node_max);
                    // printf("g:%d\n", ig);
                    graft = vector_get(tree_vec, ig);
                }
                // print_tree(root);
#ifdef DEBUG
                printf("\tPRUNE-REGRAFT p:%d (lab: %s), r:%d (lab: %s)\n", prune->id, prune->label, graft->id, graft->label);
#endif
                pr_res = prune_regraft(prune, graft, root);
                prune_res = prune;
            }
            check_subtree_losses(prune_res, tree_vec, loss_vec, k_loss, sigma, n);
//        }
    }
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
anneal(node_t *root, int sigma[], vector tree_vec, int n, int m, int k, double alpha, double beta,  int **inmatrix, double start_temp, double cooling_rate, double min_temp, double *best_lh, int MAX_LOSSES) {
    double current_temp = start_temp;

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
    for (int i = 0; i < n; i++) { current_sigma[i] = sigma[i]; }

    node_t *current_root = treecpy(root, &current_tree_vec, &current_losses_vec, current_sigma, n);

    unsigned int step = 0;
    


    double current_lh = tree_loglikelihood(current_root, tree_vec, current_sigma, inmatrix, n, m, alpha, beta);

    unsigned int not_changing_solution = 0;

    for (int phase = 1; phase <= 2; phase++) {
        if (phase == 2) {
            current_temp = 100000.0;
            cooling_rate = 0.0001;
        }
        printf("Current phase: %d\nCurrent step\t\tCurrent loglikelihood\t\tCurrent temperature\n", phase);
        while (current_temp > min_temp && not_changing_solution < 200000) {

#ifdef DEBUG
            printf("Current step %d:\n", step);
#endif
            // Create a modifiable copy
            vector copy_tree_vec;
            vector_init(&copy_tree_vec);

            vector copy_losses_vec;
            vector_init(&copy_losses_vec);
            int copy_kloss[m];
            for (int i = 0; i < m; i++) { copy_kloss[i] = current_kloss[i]; }

            int copy_sigma[n];
            for (int i = 0; i < n; i++) { copy_sigma[i] = current_sigma[i]; }

            node_t *copy_root = treecpy(current_root, &copy_tree_vec, &copy_losses_vec, copy_sigma, n);


            for (int i = 0; i < vector_total(&copy_tree_vec); i++) {
                node_t *n = vector_get(&copy_tree_vec, i);
//            printf("nodecheck: %d: lab: %s id: %d\n", i, n->label, n->id);
                assert(n->id == i);
            }
            assert(vector_total(&copy_tree_vec) == vector_total(&current_tree_vec));
            assert(vector_total(&copy_losses_vec) == vector_total(&current_losses_vec));


            neighbor(copy_root, &copy_tree_vec, copy_sigma, m, n, k, &copy_losses_vec, copy_kloss, MAX_LOSSES, phase);
//        fflush(stdout);

//        for (int i = 0; i < vector_total(&copy_tree_vec); i++) {
//            node_t *n = vector_get(&copy_tree_vec, i);
////            printf("nodecheck: %d: lab: %s id: %d\n", i, n->label, n->id);
//            if (n != NULL)
//            assert(n->id == i);
//        }
//        print_tree(copy_root);
//        print_tree(current_root);

#ifdef DEBUG
            printf("neigh\n");
#endif

//        print_tree(copy_root);

            double new_lh = tree_loglikelihood(copy_root, copy_tree_vec, copy_sigma, inmatrix, n, m, alpha, beta);

            double acceptance = accept_prob(current_lh, new_lh, current_temp);

            double p = genrand_real1();

            if (acceptance > p) {

                // Destroy current solution
                destroy_tree(current_root);

                vector_free(&current_tree_vec);
                vector_init(&current_tree_vec);

                vector_free(&current_losses_vec);
                vector_init(&current_losses_vec);


                // Copy new solution to current
                current_lh = new_lh;
                for (int i = 0; i < m; i++) { current_kloss[i] = copy_kloss[i]; }
                for (int i = 0; i < n; i++) { current_sigma[i] = copy_sigma[i]; }

#ifdef DEBUG
                printf("Accepted\n");
                for (int i = 0; i < n; i++) { printf("%02d ",current_sigma[i]); } printf("curr sigma ^ before\n");
#endif

                current_root = treecpy(copy_root, &current_tree_vec, &current_losses_vec, current_sigma, n);

#ifdef DEBUG
                for (int i = 0; i < n; i++) { printf("%02d ",current_sigma[i]); } printf("curr sigma ^ after\n");
                print_tree(current_root);
#endif

//            for (int i = 0; i < vector_total(&current_tree_vec); i++) {
//                node_t *n = vector_get(&current_tree_vec, i);
////                printf("nodecheck: %d: lab: %s id: %d\n", i, n->label, n->id);
//                assert(n->id == i);
//            }


                not_changing_solution = 0;
            } else {
                not_changing_solution++;
            }

            // Destroy neighbor
            destroy_tree(copy_root);
            vector_free(&copy_tree_vec);
            vector_free(&copy_losses_vec);


            current_temp *= (1 - cooling_rate);

            ++step;
            if (step % 1000 == 0 || step == 1) {
                printf("%d\t\t\t%lf\t\t\t%lf\n", step, current_lh, current_temp);
                // print_tree(current_root);
            }

        }
        print_tree(current_root);
    }

    *best_lh = current_lh;

    // Change SIGMA to best solution

    vector_init(&tree_vec);

    vector best_losses_vec;
    vector_init(&best_losses_vec);

    for (int i = 0; i < n; i++) { sigma[i] = current_sigma[i]; }
    node_t *best_root = treecpy(current_root, &tree_vec, &best_losses_vec, sigma, n);

    // Return best TREE
    return best_root;
}