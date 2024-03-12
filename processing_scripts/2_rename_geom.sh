#!/bin/bash

for dir in slab_atoms_*; do
    if [ -d "$dir" ]; then
        file="$dir/$dir.in"
        if [ -e "$file" ]; then
            mv "$file" "$dir/geometry.in"
            echo "Renamed $file to $dir/geometry.in"
        else
            echo "Error: $file not found in $dir"
        fi
    fi
done

