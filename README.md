Simulated Annealing Single-Cell (SASC) -- cancer progression inference
===================

Work In progress
================

Compile
--------

**OSX**
```bash
gcc -o sasc *.c -std=c99
```

**Linux**
```bash
gcc -o sasc *.c -std=gnu99 -g -lm
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

**Name of mutations file**

This optional file specifies the name of mutations (parameter `-e`). Each mutation's name must be on a different line (separated by `\n`), and the names are assigned to columns from left to right in the input file from the mutations' name file from top to bottom. If this file is not provided, mutation are progressively named from `1` to the maximum number of mutations.

An example of the mutations' name file can be seen in [pat4_mut.txt](data/scg_gawad/pat4_mut.txt).

Parameters
----------

**Required**

- `-n CELLS`: Number of cells in the input file.
- `-m MUTATIONS`: Number of mutations in the input file.
- `-k DOLLO(k)`: K value of Dollo(k) model used as phylogeny tree.
- `-a ALPHA`: False Negative rate in the input file.
- `-b BETA`: False Positive rate in the input file.
- `-i INPUT_FILE`: Path of the input file.

**Optional**
- `-d DELETIONS`: Maximum number of total deletions allowed in the solution. By default the value is set to have no restriction (+INF).
- `-e MUTATION_NAME_FILE`: Path of the input file. If this parameter is not used then the mutation will be named progressively from `1` to `MUTATIONS`

**Simulated Annealing parameters**

Although we do not recommend to change these parameters it is possible to modify them by accessing the source code.
Following is a list of the parameters used for the simulated annealing process and their position in the code.

- `START_TEMP - (line: sasc.c -> 95)`: set by default to `100000` is the starting temperature of the simulated annealing algorithm.
- `COOLING_RATE - (line: sasc.c -> 96)`: set by default  to `0.00001` is the rate of which the initial temperature decrease at each iteration.
- `MIN_TEMP - (line: sasc.c -> 97)`: set by default  to `0.0001` is the temperature at which the simulated annealing algorithm stops and return the best solution found.
- `REPETITIONS - (line: sasc.c -> 98)`: set by default  to `5` is the number of repetitions that the simulated annealing algorithm performs. The best solution is therefore the most likely from the best solutions returned by the repetitions.

The previous parameters are empirically found to be optimal, therefore we recommend to not change them but instead to try different values of `ALPHA`, `BETA` and `DOLLO(K)`.

Output
---------
Work In progress.

Run SASC
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