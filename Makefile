run: main.o mesh.o data_structures.o boundary.o solver.o vtk_export.o numerical_flux.o initial_func.o
	g++ -O3 -o run.out bin/main.o bin/solver.o bin/mesh.o bin/data_structures.o bin/boundary.o bin/vtk_export.o bin/numerical_flux.o bin/initial_func.o -fopenmp

main.o: src/main.cpp
	g++ -O3 -c -o bin/main.o src/main.cpp -fopenmp

mesh.o: src/mesh.cpp src/mesh.h
	g++ -O3 -c -o bin/mesh.o src/mesh.cpp -fopenmp

data_structures.o: src/data_structures.cpp
	g++ -O3 -c -o bin/data_structures.o src/data_structures.cpp -fopenmp

boundary.o: src/boundary.cpp
	g++ -O3 -c -o bin/boundary.o src/boundary.cpp -fopenmp

solver.o: src/solver.cpp
	g++ -O3 -c -o bin/solver.o src/solver.cpp

vtk_export.o: src/vtk_export.cpp
	g++ -O3 -c -o bin/vtk_export.o src/vtk_export.cpp

numerical_flux.o: src/numerical_flux.cpp
	g++ -O3 -c -o bin/numerical_flux.o src/numerical_flux.cpp

initial_func.o: src/initial_func.cpp src/initial_func.h
	g++ -O3 -c -o bin/initial_func.o src/initial_func.cpp

clean:
	rm bin/*.o