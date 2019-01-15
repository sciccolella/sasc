Simulated Annealing Single Cell inference (SASC) tool -- cancer progression inference
===================

SASC is  a new model and a robust framework based on Simulated Annealing for the inference of cancer progression from the SCS data.
The main objective is to overcome the limitations of the Infinite Sites Assumption by introducing a version of the k-Dollo parsimony model which indeed allows the deletion of mutations from the evolutionary history of the tumor. 

<!-- A detailed description of the framework can be found in published version of the paper [Inferring Cancer Progression from Single Cell Sequencing while allowing loss of mutations](#). -->

Compile
--------

SASC can be downloaded and compiled easily using the following commands:
```bash
git clone URL
cd FOLDER
make
```

Input Files
-------------

**Single Cell file**

The input file (specified by the `-i` parameter) is expected to be a ternary matrix file where the rows represent the cells and the columns the mutations. Each cell must be separated by a space or by a tab (`\t`). Each cell of the matrix can be:

| Value of cell | Meaning |
| ------------- | ------------- |
| I[i,j] = 0    | Mutation *j* is not observed in cell *i*  |
| I[i,j] = 1    | Mutation *j* is observed in cell *i*  |
| I[i,j] = 2    | There is no information for mutation *j* in cell *i*, i.e. low coverage  |

An example of the input file can be seen in [pat4.txt](data/scg_gawad/pat4.txt).

**Mutations file**

This optional file specifies the name of mutations (parameter `-e`). Each mutation's name must be on a different line (separated by `\n`), and the names are assigned to columns from left to right in the input file from the mutations' name file from top to bottom. If this file is not provided, mutation are progressively named from `1` to the maximum number of mutations.

An example of the mutations' name file can be seen in [pat4_mut.txt](data/scg_gawad/pat4_mut.txt).

**FN rates file**

**Prior losses file**

Usage
----------

**Input Parameters (required)**

- `-n [INT]`: Number of cells in the input file.
- `-m [INT]`: Number of mutations in the input file.
- `-k [INT]`: K value of Dollo(k) model used as phylogeny tree.
- `-a [FLOAT/STRING]`: False Negative rate in the input file or path of the file containing different FN rates for each mutations.
- `-b [FLOAT]`: False Positive rate in the input file.
- `-i [STRING]`: Path of the input file.

**Model parameters (optional)**
- `-d [INT]`: Maximum number of total deletions allowed in the solution. By default the value is set to have no restriction (+INF).
- `-e [STRING]`: Path of the input file. If this parameter is not used then the mutation will be named progressively from `1` to `MUTATIONS`
- `-g [FLOAT/STRING]`: Loss rate in the input file or path of the file containing different GAMMA rates for each mutations.
- `-r [INT]`: Set the total number of Simulated Annealing repetitions. Default is 5.

**Output parameters (optional)**
- `-l`: If this flag is used SASC will output a mutational tree with cells attached to it. Ortherwise cells will not present.
- `-x`: If this flag is used, SASC will additionally output the expected matrix E.

**Simulated Annealing parameters (optional)**
- `-S [FLOAT]`: Starting temperature of the Simulated Annealing algorithm.
- `-C [FLOAT]`: Cooling rate of the Simulated Annealing algorithm.

**Error learning  parameters (optional)**
- `-A [FLOAT]`: Standard deviation for new FN discovery.
- `-B [FLOAT]`: Standard deviation for new FP discovery.
- `-G [FLOAT]`: Standard deviation for new GAMMA discovery.

Output
---------
SASC has three different output formats that can be toggled with different arguments.

**Mutational Tree**

This is the standard output; SASC will generate a mutational tree in DOT format with no cells attached as leaves of the tree. An example of this output is shown in [Simulation mutational tree](data/results/simulation_scs_mlt.gv).

**Mutational Tree with cells as leaves**

By toggling option `-l` SASC will instead output a mutational tree in DOT format with cells attached as leaves of the tree. An example can be found at [Simulation mutational tree with cells](data/results/simulation_scs_mlt_cells.gv).

**Expected Matrix**

In addition to the previous formats SASC can output the expected matrix *E* such in [Simulation expected matrix](data/results/simulation_scs_out.txt) using the `-x` flag.

Usage examples
--------
SASC can then be run using the previously described parameters. Here we show a list of run and their results.

**Childhood Lymphoblastic Leukemia (patient 4)**
```bash
./sasc -i data/scg_gawad/pat4.txt -m 78 -n 143 -a 0.3 -b 0.001 -k 3 -d 10 -e data/scg_gawad/pat4_mut.txt 
```

The results is shown in [pat4_mlt.gv](data/results/pat4_mlt.gv), while with a post-processing procedure that collapse all the simple path it is possible to obtain a result as shown in [gawad_pat4.pdf](data/results/gawad_pat4.pdf).

**Simulation 1 of Exp. 1**
```bash
./sasc -i data/simulated/sim1_scs.txt -m 30 -n 150 -a 0.15 -b 0.001 -k 3
```

The results is shown in [sim1_scs_mlt.gv](data/results/sim1_scs_mlt.gv).