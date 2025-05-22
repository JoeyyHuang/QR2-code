#! /bin/bash
# EPW of Graphene
############################
mpirun -np 56 /opt/qe-7.3.1-qr2/bin/pw.x -npool 56 < nscf.in > nscf.out
mpirun -np 56 /opt/qe-7.3.1-qr2/bin/epw.x -npool 56 < epw.in > epw.out
