#! /bin/bash
# EPW of Graphen for Raman calculation
############################
mpirun -np 56 epw.x -npool 56 < epw-raman.in > epw-raman.out
mpirun -np 56 raman_pp.x < raman_pp1.in > raman_pp1.out
mpirun -np 56 raman_pp.x < raman_pp2.in > raman_pp2.out
