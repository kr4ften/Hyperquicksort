benchmark: main.c
	gcc -o hqs main.c -Ofast

debug: main.c
	gcc -o hqs main.c -g 

