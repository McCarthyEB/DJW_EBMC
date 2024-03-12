#!/bin/bash

#Change sbatch to qsub for PBS schedulers like that used on Young! 

for script in opt_slab_atoms_5_*.sh; do
    if [ -e "$script" ]; then
        sbatch "$script"
        #qsub "$script"
        echo "Submitted $script to queue"
    else
        echo "Error: $script not found"
    fi
done

