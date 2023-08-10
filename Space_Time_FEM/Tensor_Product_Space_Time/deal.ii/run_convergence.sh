#!/bin/bash

# for r between 0 and 3 and s between 1 and 4 (inclusive) run simulation
for r in {0..3}
do
    for s in {1..4}
    do
        # run simulation with r and s
        echo "Running simulation with r=$r and s=$s"
        ./main -r $r -s $s >> logs/output_r${r}_s${s}.txt
    done
done