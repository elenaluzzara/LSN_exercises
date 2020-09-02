CC = g++
CFLAGS = -Wall -O3 --std=c++11


2es.exe: random.o TSP_lib.o 2es.o
	mpicxx 2es.cpp random.o TSP_lib.o 
2es.o: 2es.cpp random.h TSP_lib.h
	mpicxx 2es.cpp random.o TSP_lib.o 
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TSP_lib.o: TSP_lib.cpp TSP_lib.h
	$(CC) -c TSP_lib.cpp -o TSP_lib.o $(CFLAGS)
clean:
	rm *.o *quadr* 





