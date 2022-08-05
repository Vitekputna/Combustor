#!/bin/bash

make
time ./a.out gamm_1200.mesh
paraview exp.vtk