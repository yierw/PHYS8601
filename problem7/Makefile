CC=gcc
CFLAGS= -std=c99
OBJ = thermal.o key.o ckey.o metropolis.o cmetropolis.o pb7v1.o 

run : $(OBJ)
	$(CC) -o run $(OBJ) $(CFLAGS)

clean: 
	rm -f run *.o *.dat
