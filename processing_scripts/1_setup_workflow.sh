#!/bin/bash

for geomfile in *.in; do
    string="${geomfile%.*}"
    mkdir -p "$string"
    mv "$geomfile" "$string/"
    cp opt.sh "opt_$string.sh"
    sed -i "s/test/$string/g" "opt_$string.sh"

    echo "Processing $geomfile finished!"
done

echo "All done!"

