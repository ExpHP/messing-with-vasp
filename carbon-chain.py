
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
	chain_length,    # number of links to include of the chain
	atomic_sep,      #
	layer_sep,       # separation between periodic copies on both axes
	specie = 'C',
):
	dimensions = np.array([
		chain_length * atomic_sep,
		layer_sep,
		layer_sep,
	])

	lattice = np.eye(3) * dimensions

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	for x in periodic_points(0., 1., chain_length):
		pb.add_particle(specie, [x, 1/2, 1/2], [True, True, True])

	return pb

def make_vasp_input(path, kpointdivs, comment, forbid_symmetry, **kwargs):

	pb = carbon_chain_poxcar(**kwargs)

	incar = Incar()
	incar['SYSTEM'] = comment
	incar['ALGO']   = 'Fast'
	incar['NSW']    = 1000
	incar['IBRION'] = 2
	incar['EDIFF']  = 1E-8
	if forbid_symmetry:
		incar['ISYM'] = 0

	kpoints = Kpoints.gamma_automatic(kpointdivs, [0,0,0])

	# FIXME this is all kinds of bad
	# It will include values we potentially don't want in the metadata, and prevents
	#  functions from having non-stringifiable arguments.
	# It also makes creates a dependency between the names used in the metadata file
	#  and the names of function arguments in this code
	metadata = dict(**kwargs)

	ws = VaspInput(
		poscar  = pb.poscar(comment=comment),
		potcar  = pb.potcar(functional='PBE'),
		incar   = incar,
		kpoints = kpoints,
		# additional files
		metadata = json.dumps(metadata),
	)

	ws.write_input(path)

def make_atomic_sep_study_inputs(path, seps_to_try, **kwargs):

	for trial_num, atomic_sep in enumerate(seps_to_try):
		output_dir = os.path.join(path, '{}hexagonal_scale_{}'.format(VASP_INPUT_PREFIX, trial_num))
		make_vasp_input(output_dir,
			atomic_sep   = atomic_sep,
			**kwargs
		)

	with open(os.path.join(path, SCRIPT_NAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

make_atomic_sep_study_inputs(
	path            = sys.argv[1],
	chain_length    = 1,
	seps_to_try     = np.linspace(1.200, 1.800, 12),
	forbid_symmetry = False,
	kpointdivs      = [30, 1, 1],
	layer_sep       = 15.0,
	comment         = 'Carbon Chain',
)
#kpoint_study(sys.argv[1])
