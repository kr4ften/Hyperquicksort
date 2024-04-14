benchmark: main.c
	gcc -o hqs main.c -Ofast -fopenmp

debug: main.c
	gcc -o hqs main.c -g -fopenmp

