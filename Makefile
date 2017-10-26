FLAGS=-Wall -O3

all: main.o
	@gcc -o cgSolver main.o $(FLAGS) -lm

main.o: main.c
	@gcc -o main.o -c main.c $(FLAGS)

help: 
	@echo ".cgSolver n nBanda -i <maxIter> -t <tolerancia> -o <saida>"

clean:
	@rm cgSolver *.o 

