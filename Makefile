run: main.o mesh.o data_structures.o boundary.o
	g++ -o run.out bin/main.o bin/mesh.o bin/data_structures.o bin/boundary.o -fopenmp

main.o: src/main.cpp
	g++ -c -o bin/main.o src/main.cpp -fopenmp

mesh.o: src/mesh.cpp
	g++ -c -o bin/mesh.o src/mesh.cpp -fopenmp

data_structures.o: src/data_structures.cpp
	g++ -c -o bin/data_structures.o src/data_structures.cpp -fopenmp

boundary.o: src/boundary.cpp
	g++ -c -o bin/boundary.o src/boundary.cpp -fopenmp

clean:
	rm bin/*.o