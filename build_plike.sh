#!/bin/bash
# Building PLIKE

install_plike () {
#    mkdir plike_v4
    cd plike_v4
#    wget http://www.weizmann.ac.il/home/tal/zdcf/plike_v4.0.f90
    gfortran plike_v4.0.f90 -o plike
    cd ..
}

install_plike || echo "Error installing PLIKE"
