FLAGS=-Wall -DLIKWID_PERFMON -O3 -llikwid -mavx -march=native

all: main.o
	@gcc -o cgSolver main.o -L/home/soft/likwid/lib  $(FLAGS) -lm

main.o: main.c
	@gcc -o main.o -c main.c -I/home/soft/likwid/include $(FLAGS)

help: 
	@echo ".cgSolver n nBanda -i <maxIter> -t <tolerancia> -o <saida>"

clean:
	@rm cgSolver *.o 

