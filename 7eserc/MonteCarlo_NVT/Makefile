CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) -o $@ $^  $(CFLAGS)
blkave.exe : blkave.o lib.o
	$(CC) -o $@ $^  $(CFLAGS)
autocorr.exe : autocorr.o lib.o
	$(CC) -o $@ $^  $(CFLAGS)
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)

clean :
	rm *.o Monte_Carlo_NVT.exe blkave.exe
