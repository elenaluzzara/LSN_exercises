CC = g++
CFLAGS = -Wall -O3 --std=c++11

1es.exe : 1es.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o 1es.exe seed.out averageL_circ.out bestpath_circ.out bestL_circ.out





