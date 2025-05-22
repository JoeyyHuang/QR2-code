#! /bin/bash
# Phonon of Graphene
############################
mpirun -np 56 pw.x <scf.in> scf.out
mpirun -np 56 ph.x <ph.in> ph.out
