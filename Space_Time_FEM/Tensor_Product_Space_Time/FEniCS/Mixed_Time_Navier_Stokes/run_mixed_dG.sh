#!/bin/bash

# This script runs the mixed dG method for a lot of different configurations
# and saves the results in the folder "results_mixed_dG/".

# create a counter
counter=1

for problem in "2D-3" "2D-2"
do
    echo "Running problem $problem"
    for quad in "Gauss-Lobatto" "Gauss-Legendre"
    do
        echo "  Using temporal quadrature $quad"
        for r_v in {0..2}
        do
            for r_p in {0..2}
            do

                # print counter
                echo "    $counter: Using r_v = $r_v and r_p = $r_p"

                # create folder for results with all subdirectories
                #mkdir -p results_mixed_dG/$problem/${quad}_rv_${r_v}_rp_${r_p}_dt_0.03125

                # run the code and save the output in the folder
                python3 navierstokes_mixed_dG.py --r_v $r_v --r_p $r_p --dt 0.03125 --time_dg_quadrature $quad --problem $problem #2>&1 | tee -a results_mixed_dG/$problem/${quad}_rv_${r_v}_rp_${r_p}_dt_0.03125/output.log
                ((counter++))
            done
        done
    done
done
