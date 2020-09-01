CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: 1es_p1.exe 1es_p2.exe 1es_p3.exe 

1es_p1.exe : 1es_p1.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
1es_p2.exe : 1es_p2.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
1es_p3.exe : 1es_p3.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
2es.exe : 2es.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
3es.exe : 3es.o lib.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o 1es_p1.exe 1es_p2.exe 1es_p3.exe 2es.exe 3es.exe seed.out





