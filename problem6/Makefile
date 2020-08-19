CC=gcc
CFLAGS= -std=c99
OBJ = functions.o key.o sampling.o pb6v2.o 

run : $(OBJ)
	$(CC) -o run $(OBJ) $(CFLAGS)

clean: 
	rm -f run *.o *.dat
