#!/bin/bash

for dir in slab_atoms_*; do
    if [ -d "$dir" ]; then
        # Check if the directory exists and is indeed a directory
        optrun_directories=("$dir"/optrun.*/)

        if [ ${#optrun_directories[@]} -gt 0 ]; then
            # Check if the directory exists in the source
            destination="$dir/old_geometries/"
            mkdir -p "$destination"  # Create the destination directory if it doesn't exist

            # Move each optrun directory to old_geometries
            for optrun_dir in "${optrun_directories[@]}"; do
                mv "$optrun_dir" "$destination"
                echo "Moved $optrun_dir to $destination"
            done
        else
            echo "Error: optrun.65755* directories not found in $dir"
        fi
    fi
done


