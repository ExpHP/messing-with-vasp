
from __future__ import division

import json
import math
import sys
import os

import numpy as np
import pymatgen as mg

from pymatgen.io.vaspio import Potcar, Incar, Kpoints, Poscar, VaspInput

from vaspy.poxcarbuilder import PoxcarBuilder, VelocityMode, Coords

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

# Returns n equally spaced values between start and end, centered in
#  the interval (so neither start nor end are included)
def periodic_points (start, end, n):
	res = np.linspace(start, end, n+1)

	# The interval above includes start; we want it centered
	offset = (end - start) / (2*(n))
	return res[:-1] + offset

# go go gadget "you call this unit testing?"
assert(np.linalg.norm(periodic_points(2.,3.,1) - np.array([2.5])) < 1E-10)
assert(np.linalg.norm(periodic_points(2.,3.,2) - np.array([2.25, 2.75])) < 1E-10)

def carbon_chain_poxcar (
	num_atoms,       # number of links to include of the chain
	atomicsep,       #
	layersep = None, # separation between periodic copies on both axes
	specie = 'C',
):
	if layersep is None:
		layersep = 20. * atomicsep

	lattice = np.array([
		[1.0,      0.0,      0.0],
		[0.0, layersep,      0.0],
		[0.0,      0.0, layersep],
	])

	lattice[:1] *= atomicsep

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	for x in periodic_points(0., 1., num_atoms):
		pb.add_particle(specie, [x, 1/2, 1/2], [True, False, False])

	return pb

def make_vasp_input(path, chain_length, atomic_sep, kpointdivs, comment, forbid_symmetry=True):

	pb = carbon_chain_poxcar(chain_length, atomic_sep, layersep=25.)

	incar = Incar()
	incar['SYSTEM'] = comment
	incar['ALGO']   = 'Fast'
	incar['NSW']    = 1000
	incar['IBRION'] = 2
	incar['EDIFF']  = 1E-8
	if forbid_symmetry:
		incar['ISYM'] = 0

	kpoints = Kpoints.gamma_automatic(kpointdivs, [0,0,0])

	metadata = {
		'atomic_sep': atomic_sep,
		'kpointdivs': kpointdivs,
		'chain_length': chain_length,
	}

	ws = VaspInput(
		poscar  = pb.poscar(comment=comment),
		potcar  = pb.potcar(functional='PBE'),
		incar   = incar,
		kpoints = kpoints,
		# additional files
		metadata = json.dumps(metadata),
	)

	ws.write_input(path)

def atomic_sep_study(path, chain_length, min_sep, max_sep, n):

	atomic_seps = np.linspace(min_sep, max_sep, n)

	for trial_num, atomic_sep in enumerate(atomic_seps):

		output_dir = os.path.join(path, '{}hexagonal_scale_{}'.format(VASP_INPUT_PREFIX, trial_num))

		make_vasp_input(output_dir,
			chain_length = chain_length,
			atomic_sep   = atomic_sep,
			kpointdivs   = [30, 1, 1],
			comment      = 'Carbon Chain, Length {}, sep={:0.7f}'.format(chain_length, atomic_sep),
			forbid_symmetry = True,
		)

	with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

#def kpoint_study(path):

	#atomic_sep = 1.426

	#for i in range(1, 16 + 1):
		#kpointdivs = [i, 1, 1]

		#output_dir = os.path.join(path, '{}hexagonal_kpoints_{:d}_{:d}_{:d}'.format(VASP_INPUT_PREFIX, *kpointdivs))
		#comment   = 'Carbon Hexagonal Lattice, kpointdivs={:d} {:d} {:d}'.format(*kpointdivs)

		#make_vasp_input(output_dir, atomic_sep, kpointdivs, comment)

	#with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		#f.write(SCRIPT_CONTENTS)

atomic_sep_study(sys.argv[1], 1, 1.400, 1.500, 25)
#kpoint_study(sys.argv[1])
