CC = g++
CCFLAGS = -I~/gromacs-5.1.5/spread/src/ # -ftree-vectorize


example: example.o pme-spread.o pme-malloc-free.o
	$(CC) -o $@ example.o pme-spread.o pme-malloc-free.o

pme-spread.o: pme-spread.cpp
	$(CC) -c $< -o $@

pme-malloc-free.o: pme-malloc-free.cpp
	$(CC) -c $< -o $@

example.o: example.cpp
	$(CC) -c $< -o $@
