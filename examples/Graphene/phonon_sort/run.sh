#! /bin/bash
# Phonon of Graphene
############################
mpirun -np 56 q2r.x <q2r.in> q2r.out
mpirun -np 1 phonon_sort.x < phonon_sort.in
