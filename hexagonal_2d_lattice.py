
import math

import numpy as np
import pymatgen as mg

from pymatgen.io.vaspio import Potcar, Incar, Kpoints, Poscar, VaspInput

# Turns out this is a map of  vasp_potcar_symbol ==> potcar_contents
#  e.g. 'Fe_pv' ==> (really big string)
#POTCAR_MAP = {'C': 'C'}

# TODO: As I don't think I have any control over the order in which species
#       are listed in the POSCAR, I better make certain that VASP is smart
#       enough to figure it out!

def make_hex_lattice_input(path, atomic_sep):
	symbols = ['C', 'C']
	unique_symbols = set(symbols) # for potcar

	lattice = atomic_sep * np.array([
		[3.0,            0.0,  0.],
		[1.5, math.sqrt(3)/2,  0.],
		[0.0,            0.0, 20.],
	])

	positions = [
		[0.0, 0.0, 0.0],
		[1.0, 0.0, 0.0],
	]
	dynamics = [
		[True, True, False],
		[True, True, False],
	]

	structure = mg.Structure(lattice, symbols, positions)

	potcar = Potcar(unique_symbols, functional='PBE')
	poscar = Poscar(structure, 'Hexagonal Lattice', selective_dynamics=dynamics)

	incar = Incar()
	incar['SYSTEM'] = 'Hexagonal Lattice, r={:0.3f}'.format(atomic_sep)
	incar['ALGO']   = 'F'
	incar['ENCUT']  = 250.0
	# NOTE ENCUT is in Incar.int_keys while vasp docs specify "Real"

	kpoints = Kpoints.gamma_automatic([9,9,1], [0,0,0])

	# NOTE: if not using keyword args, it seems easy to accidentally transpose the
	#       input parameters.  This error will not be detected, and it will cause
	#       the files to be written with the wrong names. :/
	ws = VaspInput(potcar=potcar, poscar=poscar, incar=incar, kpoints=kpoints)

	ws.write_input(path)

for scale in range(138,145):
	atomic_sep = scale * .01
	make_hex_lattice_input('hexagonal-{:0.3f}'.format(atomic_sep), atomic_sep)
