CC = g++
CFLAGS = -c -std=c++11 -Wall

all: prog_exec

prog_exec: main_prog.o isapprox.o ode_system_classes.o
	$(CC) main_prog.o isapprox.o ode_system_classes.o -o prog_exec

main_prog.o: main_prog.cpp calculate.h eigenvalues_calc.h integrate_system.h jacobi_am_from_boost.h linspace.h monodromy_matrix.h save_solution_function.h
	$(CC) $(CFLAGS) main_prog.cpp

isapprox.o: isapprox.cpp isapprox.h
	$(CC) $(CFLAGS) isapprox.cpp

ode_system_classes.o: ode_system_classes.cpp ode_system_classes.h
	$(CC) $(CFLAGS) ode_system_classes.cpp

clean:
	rm -rf *.o prog_exec
