#!/bin/bash

# Start simulations

# General settings
SUBCLONES=4
MUTATIONS=20
CELLS=100
#FN=0.3
#FP=0.0001
#NA=0.15
#DIR="data/simulated/subclones_$SUBCLONES-m_$MUTATIONS-n_$CELLS-fn_$FN-fp_$FP-na_$NA"
EXPERIMENT=2
DIR="data/simulated/exp$EXPERIMENT"

CORES=16

# SASC settings
SASC_FN=0.3
SASC_FP=0.01
SASC_K=2

for ID in `seq 1 50`;
do
    FILE="$DIR/sim$ID"
    FILE+="_scs.txt"
    echo $FILE

    ./sasc -i $FILE -m $MUTATIONS -n $CELLS -a $SASC_FN -b $SASC_FP -k $SASC_K &

    if [ $(($ID % $CORES)) -eq "0" ]; then
        wait
    fi
done

wait