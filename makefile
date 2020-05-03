CC=gcc
CFLAG= -Wall -pedantic -std=c99 -I.. -O3 -ffast-math -Wsuggest-attribute=const

.PHONY: all 
all: go_ai go_serial go_omp

go_ai: main_ai.c go.o
	$(CC) $(CFLAG) -o go_ai main_ai.c go.o -lm -fopenmp

go_omp: main_omp.c
	$(CC) $(CFLAG) -o go_omp main_omp.c go.o -lm -fopenmp

go_serial: main_serial.c go.o
	$(CC) $(CFLAG) -o go_serial main_serial.c go.o -lm

go.o: go.c
	$(CC) $(CFLAG) -c -o go.o go.c #-mssse3 #if you define VECTOR (also make CC=gcc44 on gpu2)


clean:
	rm -f *.o *~
