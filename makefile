kngmap: common.o main.o edlib.o
	g++ -g -fopenmp -O0 -o kngmap main.o common.o edlib.o -lz
main.o: main.cpp common.h
	g++ -g -fopenmp -c -O0 main.cpp -lz 
common.o: common.cpp common.h edlib.h
	g++ -g -fopenmp -c -O0 common.cpp
edlib.o: edlib.cpp edlib.h
	g++ -g -fopenmp -O0 -c edlib.cpp
clean:
	rm *.o kngmap
