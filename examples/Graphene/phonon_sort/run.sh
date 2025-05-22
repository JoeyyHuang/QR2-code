#! /bin/bash
# Phonon of Graphene
############################
#mpirun -np 56 /opt/qe-7.3.1-qr2/bin/q2r.x <q2r.in> q2r.out
mpirun -np 1 /opt/qe-7.3.1-qr2/bin/phonon_sort.x < phonon_sort.in
