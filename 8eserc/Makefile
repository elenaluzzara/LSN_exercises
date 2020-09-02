CC = g++
CFLAGS = -Wall -O3 --std=c++11


1es.exe : 1es.o random.o 
	$(CC) -o $@ $^  $(CFLAGS)
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o 1es.exe seed.out





