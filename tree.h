/*! \file tree.h
    \brief Include most of the operation involving Trees.
    
    Details.
*/

#include <stdbool.h>
#include "vector.h"

typedef struct Node {
  int id;
  char label[255];
  int loss;
  int mut_index;

  struct Node *first_child;
  struct Node *next_sibling;
  struct Node *previous_sibling;
  struct Node *parent;
} node_t;

/**
 * @brief Create a new node.
 * 
 * @param label Name of the node
 * @param mut_index Index of mutation encoded by the node
 * @param id Node's ID, this should be used with a global incrementing counter
 * @param tree Array of tree[ID] = *node_t, used to access nodes directly
 * @return node_t* Return pointer to new allocated node
 */
node_t * node_new(char *label, int mut_index, int id);

/**
 * @brief Append node as child of node
 * 
 * @param parent Pointer to parent
 * @param node Pointer to node
 */
void node_append(node_t *parent, node_t *node);

/**
 * @brief Deallocate a tree
 * 
 * @param node Root of the (sub)tree
 */
void destroy_tree(node_t *node);

/**
 * @brief Print tree in DOT code on stdout
 * 
 * @param root Root of the (sub)tree
 */
void print_tree(node_t *root);

/**
 * @brief Check if cand_ancestor is an ancestor of node
 * 
 * @param node Current node
 * @param cand_ancestor Candidate ancestor
 * @return true cand_ancestor is ancestor of node
 * @return false cand_ancestor is NOT ancestor of node
 */
bool is_ancestor(node_t *node, node_t *cand_ancestor);

/**
 * @brief Print tree in DOT code on stdout and attach the leaves to the nodes
 * 
 * @param root Root of the tree
 * @param tree Array of *node_t
 * @param leaves Leaves assignment (SIGMA)
 * @param MAX Number of cells
 */
void print_tree_leaves(node_t *root, node_t *tree[], int leaves[], int MAX);

/**
 * @brief Print tree in DOT code on outpath
 * 
 * @param root Root of the (sub)tree
 * @param outpath Path to output file
 */
void fprint_tree(node_t *root, char *outpath);

/**
 * @brief Print tree in DOT code on stdout and attach the leaves to the nodes
 * 
 * @param root Root of the tree
 * @param tree Array of *node_t
 * @param leaves Leaves assignment (SIGMA)
 * @param MAX Number of cells
 * @param outpath Path to output file
 */
void fprint_tree_leaves(node_t *root, vector *tree_vec, int leaves[], int MAX, char *outpath);

/**
 * @brief Get the genotype profile of the selected node
 * 
 * @param node Node to get the genotype profile
 * @param genotype Array int[MUTATIONS] initialized to all 0s
 */
void get_genotype_profile(node_t *node, int genotype[]);

/**
 * @brief Return a copy of the tree rooted in root. It also modifies the values of tree[] and max_id
 * 
 * @param root Root of the tree
 * @param tree Array of *node_t -- this will be modified accordingly
 * @param max_id Counter of assigned id -- this will be modified accordingly
 * @return node_t* Pointer to the root of the new tree
 */
node_t * treecpy(node_t *root, vector *tree_vec, vector *losses_vec, int sigma[], int n);


void node_detach(node_t *node);
bool is_loss_valid(node_t * loss);
bool is_already_lost(node_t *node, int mut_index);

void node_delete(node_t *node, vector *tree_vec, vector *loss_vec, int *k_loss, int *sigma, int n);