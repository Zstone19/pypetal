#!/bin/bash
# Building all the packages required for PETL
# The only optional arguments are u (if local install), 
#   f (the fortran compiler if specified),
#   and p (if PLIKE is installed)

cd src/petl

# Read in optional arguments
user=false
fval=false
pval=true

while getopts u:f:p: flag
do
    case "${flag}" in
        u) user=${OPTARG};;
        f) fval=${OPTARG};;
        p) pval=${OPTARG};;
    esac
done

install_plike () {
    mkdir plike_v4
    cd plike_v4
    wget http://www.weizmann.ac.il/home/tal/zdcf/plike_v4.0.f90 
    gfortran plike_v4.0.f90 -o plike
    cd .. 
}


install_javelin () {
    wget https://github.com/nye17/javelin/archive/refs/heads/master.zip
    unzip master.zip
    rm master.zip
    mv javelin-master javelin
    rm -r javelin-master

    # Compile JAVELIN
    cd javelin

    if [[ $1 = true && $2 = false ]]; then
        python setup.py install --user
    elif [[ $1 = true && $2 != false ]]; then
        python setup.py config_fc --fcompiler=$fval install --user
    elif [ $2 != false ]; then
        python setup.py config_fc --fcompiler=$fval install
    else
        python setup.py install
    fi

    cd ..
}


if [ $pval = true ]; then
    #Install PLIKE
    #NOTE: This assumes a gfortran compiler is installed
    install_plike || echo "Error installing PLIKE"
fi 

# Next, download JAVELIN
# NOTE: This requires a Fortran compiler >F90 and Git, all other packages are installed with PETL
install_javelin $user $fval || echo "Error installing JAVELIN"