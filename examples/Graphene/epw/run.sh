#! /bin/bash
# EPW of Graphene
############################
mpirun -np 56 pw.x -npool 56 < nscf.in > nscf.out
mpirun -np 56 epw.x -npool 56 < epw.in > epw.out
