#!/bin/bash 

for i in {0..20}
do
	sed '1,23d' config$i.equ > config.final #questo taglia le prime 23 righe di config$i.final perchè ci sono delle scritte e tiene solo l'orientazione dei 50 spin e la copia in config.final
	echo $i | cat - input.dat2 > input.dat #questo aggiunge la prima riga a input.dat che sarà la temperatura che però si dovrà sommare a 0.5 e dividere per 20
	./Monte_Carlo_ISING_1D.exe 
	./clean.sh #questo è per ripulire gli output.'prop'.0
done 
