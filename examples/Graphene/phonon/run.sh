#! /bin/bash
# Phonon of Graphene
############################
mpirun -np 56 /opt/qe-7.3.1-qr2/bin/pw.x <scf.in> scf.out
#mpirun -np 56 /opt/qe-7.3.1-qr2/bin/ph.x <ph.in> ph.out
