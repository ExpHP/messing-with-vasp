
from __future__ import division

import math
import sys
import os

import numpy as np
import pymatgen as mg

from pymatgen.io.vaspio import Potcar, Incar, Kpoints, Poscar, VaspInput

# Turns out this is a map of  vasp_potcar_symbol ==> potcar_contents
#  e.g. 'Fe_pv' ==> (really big string)
#POTCAR_MAP = {'C': 'C'}

# TODO: As I don't think I have any control over the order in which species
#       are listed in the POSCAR, I better make certain that VASP is smart
#       enough to figure it out!

VASP_INPUT_PREFIX = 'vi_'
SCRIPT_NAME = 'sbatch.me'
SCRIPT_CONTENTS = '''#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o out.%j
#SBATCH -J vasp-single
echo "The Job ID is $SLURM_JOB_ID"

source /cm/shared/apps/intel/bin/compilervars.sh intel64
export PATH=/cm/shared/apps/mvapich2/intel/64/1.9/bin:$PATH

for d in {}*; do
(
   cd $d
   if [ ! -e finished ]; then
     echo "Start cell optimization calculation for $d"
     srun --mpi=none vasp
     touch finished
   else
     echo "Nothing to do for $d"
   fi
)
touch finished
echo "All jobs done"
done
'''.format(VASP_INPUT_PREFIX)


def make_hex_lattice_input(path, atomic_sep):
	symbols = ['C', 'C']
	unique_symbols = set(symbols) # for potcar

	lattice = atomic_sep * np.array([
		[3.0,            0.0,  0.],
		[1.5, math.sqrt(3)/2,  0.],
		[0.0,            0.0, 20.],
	])

	positions = [
		[0.0, 0.0, 1/2],
		[1/3, 0.0, 1/2],
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

def make_test_dir(path):

	for scale in range(1380,1450):
		atomic_sep = scale * .001

		input_dir = os.path.join(path, '{}hexagonal-{:0.3f}'.format(VASP_INPUT_PREFIX, atomic_sep))

		make_hex_lattice_input(input_dir, atomic_sep)

	with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

make_test_dir(sys.argv[1])
