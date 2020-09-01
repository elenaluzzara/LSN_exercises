#!/bin/bash 

echo 0 | cat - input.gas > input.dat
./MolDyn_NVE.exe

for i in {0..6}
do
	echo 1 | cat - input.gas > input.dat
	./MolDyn_NVE.exe
done 


