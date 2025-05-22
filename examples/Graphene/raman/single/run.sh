#! /bin/bash
# EPW of Graphen for Raman calculation
############################
rm -rf *out qraman* tensor*
mpirun -np 56 /opt/qe-7.3.1-qr2/bin/epw.x -npool 56 < epw-raman.in > epw-raman.out
mpirun -np 56 /opt/qe-7.3.1-qr2/bin/raman_pp.x < raman_pp.in > raman_pp.out
