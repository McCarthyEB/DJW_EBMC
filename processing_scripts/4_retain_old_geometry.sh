#!/bin/bash

base_directory=$(pwd)

for directory in "$base_directory"/slab_atoms_5_*; do
    if [ -d "$directory" ]; then
        cd "$directory"
        if [ ! -d "old_geometries" ]; then
            mkdir old_geometries
        else
            if [ -f "geometry.in" ]; then
                mv "geometry.in" "old_geometries/geometry_original.in"
                echo "Moved geometry.in for $directory safely to old geometries! I hope!"
            fi
        fi
        cd "$base_directory"
    fi
done

