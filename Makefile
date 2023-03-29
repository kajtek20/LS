CC = g++
CFLAGS =   -std=c++1z -Wall -O2 -finline-functions   -lnfft3_threads -lfftw3_threads -lfftw3 -lm -lc -march=native -fopenmp
CFLAGS2 =  -std=c++1z -Wall -O3 -finline-functions    -lm -lc -march=native -fopenmp

# ****************************************************
# Targets needed to bring the executable up to date



main: main.o fast_LS.o fit.o interface.o
	$(CC)  -o LS_and_fit main.o fast_LS.o fit.o interface.o $(CFLAGS)

# The main.o target can be written more simply

main.o: main.cpp fast_LS.h fit.h interface.h
	$(CC) $(CFLAGS) -c main.cpp

fast_LS.o: fast_LS.cpp fast_LS.h
	$(CC) $(CFLAGS) -c fast_LS.cpp

fit.o: fit.cpp fit.h
	$(CC) $(CFLAGS2) -c fit.cpp

interface.o: interface.cpp interface.h
	$(CC) $(CFLAGS2) -c interface.cpp




.PHONY : clean
clean :
	-rm *.o LS_and_fit
