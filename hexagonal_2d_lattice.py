
from __future__ import division

import math
import sys
import os

import numpy as np
import pymatgen as mg

from pymatgen.io.vaspio import Potcar, Incar, Kpoints, Poscar, VaspInput

from vaspy.poxcarbuilder import PoxcarBuilder, VelocityMode, Coords

# Turns out this is a map of  vasp_potcar_symbol ==> potcar_contents
#  e.g. 'Fe_pv' ==> (really big string)
#POTCAR_MAP = {'C': 'C'}

# TODO: As I don't think I have any control over the order in which species
#       are listed in the POSCAR, I better make certain that VASP is smart
#       enough to figure it out!

VASP_INPUT_PREFIX = 'vc_'
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
     srun --mpi=none vasp 1>stdout 2>stderr
     touch finished
   else
     echo "Nothing to do for $d"
   fi
)
done
touch finished
echo "All jobs done"
'''.format(VASP_INPUT_PREFIX)

def hexagonal_unitcell (
	specie,
	atomicsep,
	layersep = None,
):
	if layersep is None:
		layersep = 20. * atomicsep

	lattice = np.array([
		[math.sqrt(3),    0.0,      0.0],
		[-math.sqrt(3)/2, 1.5,      0.0],
		[0.0,             0.0, layersep],
	])

	lattice[:2] *= atomicsep

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	pb.add_particle(specie, [0.0, 0.0, 0.0], [True, True, False])
	pb.add_particle(specie, [1/3, 2/3, 0.0], [True, True, False])

	return pb

# Alternate specification of a hexagonal lattice which uses an asymmetric unit cell. (a = 3, b = sqrt(3))
# Physically speaking, it is by all means identical to the standard hexagonal lattice unit cell.
#  BUT BEWARE!!!  It might catch some software off-guard...
#   
# It is here for bookeeping purposes, as it is the very first structure I tried studying in VASP.
#  When it was first used on 2015-02-15, it was discovered to produce anomalous results in VASP unless
#  symmetry optimizations were explicitly disabled (ISYM=0). (namely, the computed ground state energy
#  is greater than it should be, and has a jump discontinuity when varying atomic separation)
def hexagonal_unitcell_alt (
	specie,
	atomicsep,
	layersep = None,
):
	if layersep is None:
		layersep = 20. * atomicsep

	lattice = np.array([
		[3.0,            0.0,      0.0],
		[1.5, math.sqrt(3)/2,      0.0],
		[0.0,            0.0, layersep],
	])

	lattice[:2] *= atomicsep

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	pb.add_particle(specie, [0.0, 0.0, 1/2], [True, True, False])
	pb.add_particle(specie, [1/3, 0.0, 1/2], [True, True, False])

	return pb


# TODO delete this old mess once I test the refactored code and am convinced it still works
'''
def make_hex_lattice_input(path, atomic_sep):
	symbols = ['C', 'C']
	unique_symbols = set(symbols) # for potcar


	# lampam (old)
	lattice_lamp = np.array([
		[3.0,            0.0,  0.],
		[1.5, math.sqrt(3)/2,  0.],
		[0.0,            0.0, 20.],
	])
	# liangbo
	lattice_liang = np.array([
		[math.sqrt(3),    0.0,  0.],
		[-math.sqrt(3)/2, 1.5,  0.],
		[0.0,             0.0, 20.],
	])

	lattice = lattice_liang

	lattice[:2] = lattice[:2] * atomic_sep
	print(lattice)

	positions_lamp = [
		[0.0, 0.0, 1/2],
		[1/3, 0.0, 1/2],
	]
	positions_liang = [
		[0.0, 0.0, 1/2],
		[1/3, 2/3, 1/2],
	]

	dynamics = [
		[True, True, False],
		[True, True, False],
	]

	lattice   = lattice_lamp
	positions = positions_lamp

	structure = mg.Structure(lattice, symbols, positions)

	potcar = Potcar(unique_symbols, functional='PBE')
	poscar = Poscar(structure, 'Hexagonal Lattice', selective_dynamics=dynamics)

	incar = Incar()
	incar['SYSTEM'] = 'Hexagonal Lattice, r={:0.3f}'.format(atomic_sep)
	incar['ALGO']   = 'Fast'
	incar['NSW']    = 1000
	incar['IBRION'] = 2
	incar['EDIFF']  = 1E-8
	incar['ISYM']   = 0 # :O :O :O my worst nightmare
	# NOTE ENCUT is in Incar.int_keys while vasp docs specify "Real"

	kpoints = Kpoints.gamma_automatic([12,12,1], [0,0,0])

	# NOTE: if not using keyword args, it seems easy to accidentally transpose the
	#       input parameters.  This error will not be detected, and it will cause
	#       the files to be written with the wrong names. :/
	ws = VaspInput(potcar=potcar, poscar=poscar, incar=incar, kpoints=kpoints)

	ws.write_input(path)
'''

def make_hex_lattice_input(path, atomic_sep, kpointdivs, comment, forbid_symmetry=True):

	pb = hexagonal_unitcell('C', atomic_sep, layersep=25.)

	incar = Incar()
	incar['SYSTEM'] = comment
	incar['ALGO']   = 'Fast'
	incar['NSW']    = 1000
	incar['IBRION'] = 2
	incar['EDIFF']  = 1E-8
	if forbid_symmetry:
		incar['ISYM'] = 0

	kpoints = Kpoints.gamma_automatic(kpointdivs, [0,0,0])

	ws = VaspInput(
		poscar  = pb.poscar(comment=comment),
		potcar  = pb.potcar(functional='PBE'),
		incar   = incar,
		kpoints = kpoints,
	)

	ws.write_input(path)

def atomic_sep_study(path):

	kpointdivs = [9,9,1]

	for scale in range(1380,1450):
		atomic_sep = scale * .001

		input_dir = os.path.join(path, '{}hexagonal_scale_{:0.3f}'.format(VASP_INPUT_PREFIX, atomic_sep))
		comment   = 'Carbon Hexagonal Lattice, r={:0.3f}'.format(atomic_sep)

		make_hex_lattice_input(input_dir, atomic_sep, kpointdivs, comment)

	with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

def kpoint_study(path):

	atomic_sep = 1.42

	for i in range(1, 16 + 1):
		kpointdivs = [i, i, 1]

		input_dir = os.path.join(path, '{}hexagonal_kpoints_{:d}_{:d}_{:d}'.format(VASP_INPUT_PREFIX, *kpointdivs))
		comment   = 'Carbon Hexagonal Lattice, kpointdivs={:d} {:d} {:d}'.format(*kpointdivs)

		make_hex_lattice_input(input_dir, atomic_sep, kpointdivs, comment)

	with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		f.write(SCRIPT_CONTENTS)


atomic_sep_study(sys.argv[1])
#kpoint_study(sys.argv[1])
