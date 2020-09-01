#!/bin/bash 

echo 0 | cat - input.solid > input.dat
./MolDyn_NVE.exe

for i in {0..8}
do
	echo 1 | cat - input.solid > input.dat
	./MolDyn_NVE.exe
done 
