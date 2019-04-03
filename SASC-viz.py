#     MIT License
#
#     Copyright (c) 2017-2019 Simone Ciccolella
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy
#     of this software and associated documentation files (the "Software"), to deal
#     in the Software without restriction, including without limitation the rights
#     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#     copies of the Software, and to permit persons to whom the Software is
#     furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in all
#     copies or substantial portions of the Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#     SOFTWARE.
#

import argparse
from colour import Color
import re
import sys

def load_from_sasc(filepath, cell_labels, show_support=False, score=True, show_color=True):
    tree = []
    with open(filepath, 'r') as f:
        tree = f.readlines()

    tot_cells = 0

    ROOT = Node(0)
    ROOT.tot_cells = len(cell_labels) + 1
    TREE = Tree(ROOT)
    if score:
        gv = tree[1:-3]
    else:
        gv = tree[1:-1]

    for line in gv:
        line = line.strip()[:-1]
        if ' -> ' in line:
            s, e = line.split(' -> ')
            
            s = s.replace('"', '')
            e = e.replace('"', '')

            if e in cell_labels:
                # this is a cell name
                s = int(s)
                x = TREE.get_node(s)
                x.support += 1
                continue
            
            s = int(s)
            try:
                e = int(e)
            except:
                sys.exit('The node "{}" is probably a cell, but it was not found in the labels provided. Please check the input'.format(e))

            x = TREE.get_node(s)
            if not x:
                x = Node(s)
                x.tot_cells = len(cell_labels) + 1
                x.show_support = show_support
                x.show_color = show_color
                TREE.add_node(x)

            y = TREE.get_node(e)
            if not y:
                y = Node(e)
                y.tot_cells = len(cell_labels) + 1
                y.show_support = show_support
                y.show_color = show_color
                TREE.add_node(y)

            TREE.add_edge(s, e)

        if 'label' in line:
            id, label = line.split(' [')

            id = int(id.replace('"', '').strip())
            label = label[: -1]
            if 'color=indianred1' in label:
                m = re.search(r'label="(.+)"', label)
                mut = str(m.group(1))
                TREE.set_deletion(id)
                TREE.set_mutations(id, [mut])
            else:
                if id != 0:
                    m = re.search(r'label="(.+)"', label)
                    mut = str(m.group(1))
                else:
                    mut = 'germline'
                TREE.set_mutations(id, [mut])

    return TREE

class Node:
    def __init__(self, id):
        self.id = id
        self.mutations = None
        self.deletion = False

        self.parent = None
        self.children = []

        self.depth = 0

        self.c_grad = None

        self.support = 0
        self.cumulative_support = 0
        self.downstream_support = 0

        self.tot_cells = 0
        self.show_support = False
        self.show_color = False
    
    def get_name(self, sep=','):
        if not self.deletion:
            return sep.join(self.mutations)
        else:
            return sep.join('%s-' % x for x in self.mutations)

    def get_s(self):
        return int((self.downstream_support / (self.tot_cells - self.parent.cumulative_support))*100)

    def print_node_dot(self, sep=','):
        c_red = Color("#FF1919")
        c_gradient = list(c_red.range_to(Color("#397D02"), 100))

        if self.parent:
            print('\t"%s" -> "%s";' % (self.parent.id, self.id))

            # s = int((self.downstream_support / (self.tot_cells - self.parent.cumulative_support))*100)
            s = self.get_s()
            
            if self.parent.parent:
                parent_s = self.parent.c_grad
                increment = int((100 - parent_s) * (s / 100.0))
                color = c_gradient[parent_s + increment]
                self.c_grad = parent_s + increment
            else:
                # print(s)
                color = c_gradient[s]
                self.c_grad = s

            print_label = None
            if self.show_support:
                print_label = '{} [s={}%]'.format(
                    self.get_name(sep=sep),
                    s
                )
            else:
                print_label = self.get_name(sep=sep)

            if self.show_color:
                color_label = ', color="{}"'.format(color)
            else:
                color_label = ''

            if not self.deletion:
                print('\t"{0}" [label="{1}"{2}];'.format(
                    self.id, 
                    print_label,
                    color_label
                    ))
            else:
                print('\t"{0}" [label="{1}"{2}, fillcolor=indianred1, style=filled];'.format(
                    self.id, 
                    print_label,
                    color_label
                    ))
    
    def calc_cumalitve_sup(self):
        if self.parent:
            self.cumulative_support = self.parent.cumulative_support + self.support
        else:
            self.cumulative_support = self.support

class Tree:
    def __init__(self, root):
        self.root = root
        self.id_to_node = {}
        self.mut_to_node = {}
        self.deletions = []
        self.edges = []

        self.add_node(root)
    
    def add_node(self, node):
        self.id_to_node[node.id] = node

    def remove_node(self, node):
        self.id_to_node.pop(node.id)
        for m in node.mutations:
            self.mut_to_node.pop(m)
        
        node.parent.children.remove(node)
        for child in node.children:
            child.parent = node.parent
            node.parent.children.append(child)
    
    def pop_node(self, node):
        self.id_to_node.pop(node.id)
        for m in node.mutations:
            self.mut_to_node.pop(m)
            
        node.parent.children.remove(node)
        for child in node.children:
            child.parent = node.parent

    def merge_nodes(self, merged, to_merge):
        if not to_merge.parent == merged:
            sys.exit('Merge two non-parent-child nodes')

        for m in to_merge.mutations:
            merged.mutations.append(m)
        self.remove_node(to_merge)
        for m in merged.mutations:
            self.mut_to_node[m] = merged
        
        merged.support += to_merge.support


    def get_node(self, id):
        try:
            return self.id_to_node[id]
        except:
            return None

    def add_edge(self, start_id, end_id):
        s_node = self.get_node(start_id)
        e_node = self.get_node(end_id)

        e_node.parent = s_node
        s_node.children.append(e_node)
        self.edges.append((start_id, end_id))
        e_node.cumulative_support = s_node.cumulative_support

    def set_mutations(self, id, mutations):
        self.get_node(id).mutations = mutations
        for mut in mutations:
            if not self.get_node(id).deletion:
                self.mut_to_node[mut] = self.get_node(id)

    def set_deletion(self, id):
        self.get_node(id).deletion = True
        self.deletions.append(self.get_node(id))

    def is_ancestor(self, anc_mut, node_mut):
        anc = self.mut_to_node[anc_mut]
        node = self.mut_to_node[node_mut]
        
        p = node.parent
        while p:
            if p == anc:
                return True
            
            p = p.parent
        return False

    def get_deletions_name(self):
        ret = []
        for d in self.deletions:
            ret.append(d.get_name())
        return ret

def collapse_simple_paths(tree, node):
    if len(node.children) == 0:
        return
    else:
        if not node.deletion:
            if len(node.children) == 1 and node.parent and not node.children[0].deletion:
                # collapse parent with its child
                only_child = node.children[0]

                tree.merge_nodes(node, only_child)
                collapse_simple_paths(tree, node)
        for child in node.children:
            collapse_simple_paths(tree, child)

def collapse_low_support(tree, node, support_th):
    if node.parent and not node.deletion and not node.parent.deletion:
        if node.get_s() < support_th and not node.parent == tree.root:
            # collapse node with its parent
            tree.merge_nodes(node.parent, node)
            
            node = node.parent
    for child in node.children:
        collapse_low_support(tree, child, support_th)

def delete_subtree(tree, node):
    if len(node.children) == 0:
        tree.pop_node(node)
    else:
        for child in node.children:
            delete_subtree(tree, child)
        tree.pop_node(node)

def __print_tree(node, ds_filter=0, sep=','):
    if len(node.children) == 0:
        if node.cumulative_support >= ds_filter:
            node.print_node_dot(sep=sep)
    else:
        if node.cumulative_support >= ds_filter:
            node.print_node_dot(sep=sep)
        for child in node.children:
            __print_tree(child, ds_filter=ds_filter, sep=sep)

def print_dot_tree(node, ds_filter=0, sep=',', show_support=False):
    print('digraph phylogeny {')
    print('\tnode [penwidth=2];')
    if show_support:
        support = ' [{} cells]'.format(node.support)
    else:
        support = ''
    print('\t"{0}" [label="{1}{2}"];'.format(node.id, ','.join(node.mutations), support))
    __print_tree(node, ds_filter=ds_filter, sep=sep)
    print('}')

def calc_supports(node, level_count):
    if node.parent:
        node.depth = node.parent.depth + 1
    level_count[node.depth] += 1
    if len(node.children) == 0:
        node.calc_cumalitve_sup()
        node.downstream_support = node.support
    else:
        node.calc_cumalitve_sup()
        for child in node.children:
            calc_supports(child, level_count)
        node.downstream_support = sum([x.downstream_support for x in node.children]) + node.support





parser = argparse.ArgumentParser(description='SASC visualitation tool', add_help=True)

parser.add_argument('-t', '--tree', action='store', type=str, required=True,
                    help='path of the input file.')
group = parser.add_mutually_exclusive_group()
group.add_argument('-E', '--cellnames', action='store', type=str,
                    required=False,
                    help="path to the cell labels file")
group.add_argument('-n', '--totcell', action='store', type=int,
                    required=False,
                    help="total number of cells")

parser.add_argument('--show-support', action='store_true', required=False,
                    help="Show the support for each node.")

parser.add_argument('--show-color', action='store_true', required=False,
                    help="Enable coloring of nodes.")

parser.add_argument('--collapse-support', action='store', type=float,
                    required=False,
                    help="Collapse path with lower support")

parser.add_argument('--collapse-simple', action='store_true',
                    required=False,
                    help="Collapse simple paths")

parser.add_argument('--sep', action='store', default=',', type=str,
                    help="Labels' separator")

args = parser.parse_args()

cells_labels = set()

if args.cellnames:
    with open(args.cellnames) as fin:
        for line in fin:
            cells_labels.add(line.strip())
else:
    for x in range(args.totcell):
        cells_labels.add('cell{}'.format(x+1))

x = load_from_sasc(args.tree, cells_labels, show_support=args.show_support, show_color=args.show_color, score=True)

from collections import defaultdict

lev_count = defaultdict(int)
lev_count[0] = 1

calc_supports(x.root, lev_count)

if args.collapse_simple:
    collapse_simple_paths(x, x.root)

if args.collapse_support:
    collapse_low_support(x, x.root, args.collapse_support)


lev_count = defaultdict(int)
lev_count[0] = 1

calc_supports(x.root, lev_count)

print_dot_tree(x.root, sep=args.sep, show_support=args.show_support)