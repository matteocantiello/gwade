#!/bin/bash

mesa-21

echo "Compiling..."
cd template
./mk &> ../mk.out
cd ..
echo "Done compiling."

echo "Starting runs."
mkdir runs

for j in 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8

do
    rm -r runs/$j
    cp -R template runs/$j    
    cd runs/$j
    sed -i "s/MMM/$j/" inlist_project
    echo "Running Model..." $j
    ./rn &> rn.out
    cd ../../
done
