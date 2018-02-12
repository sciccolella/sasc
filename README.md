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

Input File(s)
-------------
Work In Progress

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
- `-e MUTATION_NAME_FILE`: Path of the input file. If this parameter is not used then the mutation will be named progressively from `1` to `MUTATIONS`

**Simulated Annealing parameters**
Although we do not recommend to change these parameters it is possible to modify them by accessing the source code.
Following is a list of the parameters used for the simulated annealing process and their position in the code.

- `START_TEMP - (line: XX)`: set initially to `100000` is the starting temperature of the simulated annealing algorithm
- `COOLING_RATE - (line: XX)`: set initially to `0.00001` is the rate of which the initial temperature decrease at each iteration


Run
--------
Work In Progress

```bash
    ./sasc -i data/simulated/subclones_7-m_20-n_100-fn_0.3-fp_0.0001-na_0.15/sim2_scs.txt -m 20 -n 100 -a 0.3 -b 0.001 -k 2

    ./sasc -i data/gawad2.txt -m 16 -n 115 -a 0.3 -b 0.1 -k 2 -e data/gawad2_mut.txt
    ./sasc -i data/hou.txt -m 18 -n 58 -a 0.43 -b 0.0000604 -k 2 -e data/hou_mut.txt
```