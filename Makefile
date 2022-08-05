run: main.o mesh.o data_structures.o boundary.o
	g++ main.o mesh.o data_structures.o boundary.o

main.o: src/main.cpp src/*.h
	g++ -c src/main.cpp

mesh.o: src/mesh.cpp src/*.h
	g++ -c src/mesh.cpp

data_structures.o: src/data_structures.cpp src/*.h
	g++ -c src/data_structures.cpp

boundary.o: src/boundary.cpp src/*.h
	g++ -c src/boundary.cpp

clean:
	rm *.o