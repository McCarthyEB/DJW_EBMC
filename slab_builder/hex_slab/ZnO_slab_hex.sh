#!/bin/bash -i

module purge
#module load python
module load python/3.7.7-intel2020u1

python3 ZnO_slab_hex.py
