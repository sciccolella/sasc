#include "tree.h"
#include "sastep.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "vector.h"

node_t *
node_new(char *label, int mut_index, int id) {
    node_t *node = malloc(sizeof(node_t));
    node->id = id;
    node->mut_index = mut_index;
    
    strcpy(node->label, label);
    node->loss = 0;
    node->first_child = NULL;
    node->next_sibling = NULL;
    node->previous_sibling = NULL;
    node->parent = NULL;

    return node;
}

void
node_append(node_t *parent, node_t *node){
    if (parent->first_child != NULL)
        parent->first_child->previous_sibling = node;
    node->next_sibling = parent->first_child;
    parent->first_child = node;
    node->parent = parent;
}

bool
is_ancestor(node_t *node, node_t *cand_ancestor) {
    node_t *par = node->parent;

    while (par != NULL) {
        if (par == cand_ancestor)
            return true;
        par = par->parent;
    }
    return false;
}

bool
is_already_lost(node_t *node, int mut_index) {
    node_t *par = node->parent;

    while (par != NULL) {
        if (par->loss == 1 && par->mut_index == mut_index)
            return true;
        par = par->parent;
    }
    return false;
}

bool
is_loss_valid(node_t * loss) {
    node_t *par = loss->parent;

    while (par != NULL) {
        if (par->mut_index == loss->mut_index)
            return true;
        par = par->parent;
    }
    return false;
}

void
node_detach(node_t *node) {
    if (node->next_sibling == NULL) {
        if (node->previous_sibling != NULL){
            node->previous_sibling->next_sibling = NULL;
            node->previous_sibling = NULL;
        } else
            node->parent->first_child = NULL;
    } else if (node->previous_sibling == NULL) {
        if (node->next_sibling != NULL) {
            node->parent->first_child = node->next_sibling;
            node->next_sibling->previous_sibling = NULL;
            node->next_sibling = NULL;
        } else
            node->parent->first_child = NULL;
    } else {
        node->previous_sibling->next_sibling = node->next_sibling;
        node->next_sibling->previous_sibling = node->previous_sibling;
        node->previous_sibling = NULL;
        node->next_sibling = NULL;
    }
}

void
node_delete(node_t *node, vector * tree_vec, vector *loss_vec, int *k_loss, int *sigma, int n) {
    
    node_t *par = node->parent;

    node_detach(node);

    int del_id = node->id;

    assert(vector_total(tree_vec) > 0);
    assert(vector_total(loss_vec) > 0);
    assert(vector_get(tree_vec, del_id) == node);

    vector_set(tree_vec, del_id, NULL);

    // If there was a cell assigned to the deleted node then re-assign it.
    for (int i =0; i < n; i++) {
        if (sigma[i] == del_id) {
            int x = random_assignment(vector_total(tree_vec) -1);
            node_t *new_node = NULL;
            while (new_node == NULL) {
                x = random_assignment(vector_total(tree_vec) -1);
                new_node = vector_get(tree_vec, x);
            }
            sigma[i] = x;
        }
    }

    k_loss[node->mut_index] -= 1;

    // Delete from loss vector
    for (int i = 0; i < vector_total(loss_vec); i ++) {
        node_t *x = vector_get(loss_vec, i);
        if (x == node) {
            vector_delete(loss_vec, i);
            break;
        }
    }

    node_t *child = node->first_child;

    while (child != NULL) {
        node_t *ns = child->next_sibling;

        node_detach(child);
        node_append(par, child);

        child = ns;
    }

    free(node);
}

void
destroy_tree(node_t *node){
    if (node == NULL) return;
    destroy_tree(node->first_child);
    destroy_tree(node->next_sibling);
    free(node);
}

 void
 print_tree_rec(node_t *node) {
     if (node == NULL) return;
     print_tree_rec(node->first_child);
     print_tree_rec(node->next_sibling);
     if(node->parent != NULL)
         printf("\t\"%d\" -> \"%d\";\n",
                 node->parent->id,
                 node->id);

     if (node->loss == 1)
          printf("\t\"%d\" [color=indianred1, style=filled, label=\"%s\"];\n",
                 node->id, node->label);
     else
         printf("\t\"%d\" [label=\"%s\"];\n",
                 node->id, node->label);
 }

 void
 print_tree(node_t *root, double score) {
   printf("%s\n", "digraph g {");
   print_tree_rec(root);
   printf("\tlabelloc=\"t\";\n\tlabel=\"Confidence score: %lf\";\n}\n", score);
 }

// void
// print_tree_leaves(node_t *root, node_t *tree[], int leaves[], int MAX) {
//     printf("%s\n", "digraph g {");
//     print_tree_rec(root);

//     for (int i = 0; i < MAX; i++) {
//         printf("\t\"%s\" -> cell%d;\n",
//                     tree[leaves[i]]->label, i+1);
//         printf("\tcell%d [shape=box]\n", i+1);
//     }

//     printf("}\n");

// }

void
fprint_tree_rec(node_t *node, FILE *fo) {
    if (node == NULL) return;
    fprint_tree_rec(node->first_child, fo);
    fprint_tree_rec(node->next_sibling, fo);

    if(node->parent != NULL)
        fprintf(fo, "\t\"%d\" -> \"%d\";\n",
               node->parent->id,
               node->id);

    if (node->loss == 1)
        fprintf(fo, "\t\"%d\" [color=indianred1, style=filled, label=\"%s\"];\n",
               node->id, node->label);
    else
        fprintf(fo, "\t\"%d\" [label=\"%s\"];\n",
               node->id, node->label);
}

void
fprint_tree(node_t *root, char *outpath, double score) {
    FILE *fp;
    fp = fopen(outpath, "w+");
    
    fprintf(fp, "%s\n", "digraph g {");
    fprint_tree_rec(root, fp);
    fprintf(fp, "\tlabelloc=\"t\";\n\tlabel=\"Confidence score: %lf\";\n}\n", score);
    fclose(fp);
}

void
fprint_tree_leaves(node_t *root, vector *tree_vec, int sigma[], int MAX, char *outpath, double score, char cell_names[][50]) {
    FILE *fp;
    fp = fopen(outpath, "w+");

    fprintf(fp, "%s\n", "digraph g {");
    fprint_tree_rec(root, fp);

    for (int i = 0; i < MAX; i++) {
        node_t *node = vector_get(tree_vec, sigma[i]);
        fprintf(fp, "\t\"%d\" -> \"%s\";\n",
                    node->id, cell_names[i]);
        fprintf(fp, "\t\"%s\" [shape=box];\n", cell_names[i]);
    }

    fprintf(fp, "\tlabelloc=\"t\";\n\tlabel=\"Confidence score: %lf\";\n}\n", score);
    fclose(fp);

}

void get_genotype_profile(node_t *node, int genotype[]) {
    if (node->mut_index == -1)
        return;
    if (node->loss == 0)
        genotype[node->mut_index] += 1;
    else
        genotype[node->mut_index] -= 1;

    get_genotype_profile(node->parent, genotype);
       
}

node_t *
nodecpy(char *label, int mut_index, int id, int loss) {
    node_t *node = malloc(sizeof(node_t));
    node->id = id;
    node->mut_index = mut_index;
    
    strcpy(node->label, label);
    node->loss = loss;
    node->first_child = NULL;
    node->next_sibling = NULL;
    node->previous_sibling = NULL;
    node->parent = NULL;

    return node;
}

void
rec_treecpy(node_t *node, node_t *cnode, vector *tree_vec, vector *losses_vec, int n, bool changed_mask[]) {
    node_t *child = node->first_child;
    if (child == NULL)
        return;

    while (child != NULL) {
        int new_id = vector_total(tree_vec);

        node_t *copy = nodecpy(child->label, child->mut_index, new_id, child->loss);
        vector_add(tree_vec, copy);

        if (copy->loss == 1)
            vector_add(losses_vec, copy);

        // if (new_id != child->id) {
        //     for (int i =0; i < n; i++) {
        //         if (sigma[i] == child->id && !changed_mask[i]) {
        //             sigma[i] = new_id;
        //             changed_mask[i] = true;
        //         } else {
        //         }
        //     }
        // }

        node_t *ncheck = vector_get(tree_vec, new_id);
        assert(ncheck->id == new_id);
        assert(ncheck == copy);

        assert(child->mut_index == copy->mut_index);

        node_append(cnode, copy);
        rec_treecpy(child, copy, tree_vec, losses_vec, n, changed_mask);

        child = child->next_sibling;
    }
}

node_t *
treecpy(node_t *root, vector *tree_vec, vector *losses_vec, int n) {
    if (vector_total(tree_vec) != 0) {
        fprintf(stderr, "ERROR: tree_vec incorrectly initialized in TREECPY.\n");
        exit(EXIT_FAILURE);
    }
    if (vector_total(losses_vec) != 0) {
        fprintf(stderr, "ERROR: losses_vec incorrectly initialized in TREECPY.\n");
        exit(EXIT_FAILURE);
    }

    node_t *croot = nodecpy(root->label, root->mut_index, 0, root->loss);

    assert(croot->id == 0);
    assert(croot->mut_index == -1);

    bool changed_mask[n];
    for (int i = 0; i < n; i++) { changed_mask[i] = false;}

    vector_add(tree_vec, croot);
    node_t *child = root->first_child;

    while (child != NULL) {
        int new_id = vector_total(tree_vec);

        node_t *copy = nodecpy(child->label, child->mut_index, new_id, child->loss);
        vector_add(tree_vec, copy);

        node_t *ncheck = vector_get(tree_vec, new_id);
        assert(ncheck->id == new_id);
        assert(ncheck == copy);


        if (copy->loss == 1)
            vector_add(losses_vec, copy);

        // if (new_id != child->id) {
        //     for (int i =0; i < n; i++) {
        //         if (sigma[i] == child->id && !changed_mask[i]) {
        //             sigma[i] = new_id;
        //             changed_mask[i] = true;
        //         } else {
        //         }

        //     }
        // }

        assert(child->mut_index == copy->mut_index);

        node_append(croot, copy);
        rec_treecpy(child, copy, tree_vec, losses_vec, n, changed_mask);

        child = child->next_sibling;
    }

    return croot;
}
