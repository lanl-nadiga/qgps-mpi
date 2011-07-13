CFLAGS=-std=c99 -O2
LD_FLAGS=-lfftw3_mpi -lfftw3 -lmpi -lm

SOURCE=$(wildcard *.c)
OBJECTS=$(SOURCE:.c=.o)

ogps-mpi: $(OBJECTS)
	$(CC) $(OBJECTS) $(LD_FLAGS) -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o ogps-mpi
