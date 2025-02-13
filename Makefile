CC=gcc
CFLAGS=-O9
LDFLAGS=
EXEC=main

all: $(EXEC)

main: partition.o main.o
	$(CC) -o recpart partition.o main.o $(CFLAGS)

bisection.o: partition.c
	$(CC) -o partition.o -c partition.c $(CFLAGS)

main.o: main.c
	$(CC) -o main.o -c main.c $(CFLAGS)

clean:
	rm *.o recpart

