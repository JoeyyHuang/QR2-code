#! /bin/bash
# EPW of Graphen for Raman calculation
############################
mpirun -np 56 epw.x -npool 56 < epw-raman.in > epw-raman.out
mpirun -np 56 raman_pp.x < raman_pp.in > raman_pp.out
