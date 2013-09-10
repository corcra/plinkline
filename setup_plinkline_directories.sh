#!/bin/bash

echo "Creating a load of folders!"

mkdir ./idlists
mkdir ./ihshap_files
mkdir ./ihsmap_files
mkdir ./test_statistics
mkdir ./test_statistics/xpehh
mkdir ./test_statistics/ihs
mkdir ./regions

for i in {1..39}
do
    mkdir ./ihshap_files/chr$i/
    mkdir ./ihmap_files/chr$i/
    mkdir ./test_statistics/xpehh/chr$i/
    mkdir ./test_statistics/ihs/chr$i/
done
