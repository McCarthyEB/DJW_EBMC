#!/bin/bash

#This is only useful and necessary for Hawk, as on Young and Archer2 the jobs are run from the scratch and work directories!

base_directory=$(pwd)

for slab_directory in "$base_directory"/slab_atoms_5_*; do
    if [ -d "$slab_directory" ]; then
        cd "$slab_directory"
        number=$(basename "$slab_directory" | sed -n 's/.*slab_atoms_5_\([0-9]*\)/\1/p')
        optrun_directory=$(find . -maxdepth 1 -type d -name "optrun.*" -print -quit)
        if [ -n "$optrun_directory" ]; then
            cd "$optrun_directory"
            job_number=$(basename "$PWD" | sed -n 's/optrun\.\([0-9]*\)/\1/p')
            source_directory="/scratch/c.c1809139/$job_number/slab_atoms_5_$number/optrun.$job_number/"
            if [ -d "$source_directory" ]; then
                cp -r "$source_directory"* .
            else
                echo "Source directory not found: $source_directory"
            fi
            cp "geometry.in" ..
            cp "qn.pckl" ..
            cd ..
        else
            echo "optrun.<jobnumber> directory not found in $slab_directory"
        fi
        cd "$base_directory"
    fi
done

