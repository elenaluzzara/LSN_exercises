#!/bin/bash 

for i in {0..20}
do

	echo $i | cat - input.dat2 > input.dat #questo aggiunge la prima riga a input.dat che sarà la temperatura che però si dovrà sommare a 0.5 e dividere per 20
	./Monte_Carlo_ISING_1D.exe > config$i.equ #questo serve per equilibrare e salvare in config$i.final la configurazione finale

done 

rm -rf *temp*
rm -rf output.*
