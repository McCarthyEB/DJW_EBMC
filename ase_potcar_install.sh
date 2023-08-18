#!/bin/bash -i
#
# ASE builds the POTCAR file itself from the set available from this path:
#               The pseudopotentials are expected to be in:
#               LDA:  $VASP_PP_PATH/potpaw/
#               PBE:  $VASP_PP_PATH/potpaw_PBE/
#               PW91: $VASP_PP_PATH/potpaw_GGA/
#
export VASP_PP_PATH=$HOME/progs/vasp

# make the required top directory
#
mkdir -p $VASP_PP_PATH
#
# move the potcar tar file
gunzip potcar_PBE.tar.gz
mv potcar_PBE.tar $VASP_PP_PATH
#
# go there and unpack     
cd  $VASP_PP_PATH
tar -xvf potcar_PBE.tar
#
#remove tar file
rm potcar_PBE.tar
