#!/usr/bin/env python3
# python 2 sucks and you know it

# 2015 Apr 8
# Hexagonal cell gen, for use in VASP

import math
import sys
from workspace_maker import VaspWorkspaceMaker
from workspace_maker import write_composite_potcar

# TODO figure out what this actually is
CARBON_SEP_ANGSTROMS = 1.42

# TODO POTCAR files must be gathered from vasp

class MyWorkspaceMaker(VaspWorkspaceMaker):
	# atomic_separation:  Distance between any two adjacent atoms in the lattice
	def __init__ (self, title, atomic_separation, dynamics = (False, False, False)):
		self._title = title
		self._dynamics = dynamics
		self._atomic_separation = float(atomic_separation)

	def write_poscar (self, f):
		# Title
		f.write("{}\n".format(self._title))

		# global scale factor
		f.write(" {:15f}\n".format(self._atomic_separation))

		# lattice cell vectors
		self._write_cubic_lattice_vectors(f, 3., math.sqrt(3.), 20.)

		particles = [
			(1./6., 0.0, 0.5),
			(2./6., 0.5, 0.5),
			(4./6., 0.5, 0.5),
			(5./6., 0.0, 0.5),
		]

		# number of atoms
		f.write(' {:d}\n'.format(len(particles)))

		f.write('Selective dynamics\n')
		f.write('Direct\n')

		# positions of atoms
		for pos in particles:
			self._write_position_with_dynamics(f, pos, self._dynamics)

	def _write_cubic_lattice_vectors (self, f, xdim, ydim, zdim):
		self._write_lattice_vector(f, xdim,  0.0,  0.0)
		self._write_lattice_vector(f,  0.0, ydim,  0.0)
		self._write_lattice_vector(f,  0.0,  0.0, zdim)

	def _write_lattice_vector (self, f, x, y, z):
		f.write((" {:.15f}"*3 + '\n').format(x,y,z))

	def _write_position_with_dynamics(self, f, pos, dynamics):
		assert len(pos) == len(dynamics) == 3
		pos_strs = ['{:.15f}'.format(x) for x in pos]
		dynamics_strs = [('T' if x else 'F') for x in dynamics]
		full_str = ' '.join(pos_strs + dynamics_strs) + '\n'
		f.write(full_str)

	def write_incar (self, f):
		f.write('SYSTEM = {}\n'.format(self._title))
		f.write('ISMEAR = 0\n') # gaussian smearing (of partial occupancies)

	def write_potcar (self, f):
		write_composite_potcar(f, ['C'])

	def write_kpoints (self, f):
		f.write("hi i'm kpoints\n")
		f.write(' 0\n') # 0 is automatic generation
		f.write('Gamma point\n') # monkhurst or gamma point
		f.write('9 9 1\n') # mesh parameter (# intersections along each axis)
		f.write('0 0 0\n') # shift

def main():
	maker = MyWorkspaceMaker("Single Hexagonal Cell of Something", CARBON_SEP_ANGSTROMS)
	maker.make_workspace(sys.argv[1])

if __name__ == '__main__':
	main()
