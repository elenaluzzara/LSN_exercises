CC = g++
CFLAGS = -Wall -O3 --std=c++11

1es.exe : 1es.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
2es.exe : 2es.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o 1es.exe 2es.exe seed.out





