
from __future__ import division

import json
import math
import sys
import os

import numpy as np
import pymatgen as mg

from pymatgen.io.vaspio import Potcar, Incar, Kpoints, Poscar, VaspInput

from vaspy.poxcarbuilder import PoxcarBuilder, VelocityMode, Coords

import argparse

TRIAL_DIR_LIST_FILENAME = 'trial_dirs'
SCRIPT_FILENAME = 'sbatch.me'
SCRIPT_CONTENTS = '''#!/bin/bash
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -o out.%j
#SBATCH -J vasp-single
echo "The Job ID is $SLURM_JOB_ID"

source /cm/shared/apps/intel/bin/compilervars.sh intel64
export PATH=/cm/shared/apps/mvapich2/intel/64/1.9/bin:$PATH

for d in `cat "{}"`; do
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
'''.format(TRIAL_DIR_LIST_FILENAME)

# Returns n equally spaced values between start and end, centered in
#  the interval (so neither start nor end are included)
def periodic_points(start, end, n):
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
	], dtype=float)

	lattice = np.eye(3) * dimensions

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	for x in periodic_points(0., 1., chain_length):
		pb.add_particle(specie, [x, 1/2, 1/2], [True, True, True])

	return pb

def make_vasp_input(path, *, kpointdivs, comment, forbid_symmetry, functional, perturb_dist, **kwargs):

	pb = carbon_chain_poxcar(**kwargs)

	incar = Incar()
	incar['SYSTEM'] = comment
	incar['ALGO']   = 'Fast'
	incar['NSW']    = 1000
	incar['NPAR']   = 4
	incar['IBRION'] = 2
	incar['EDIFF']  = 1E-8
	incar['EDIFFG'] = -0.001
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
		poscar = pb.poscar(
			perturb_dist=perturb_dist,
			comment=comment,
		),
		potcar = pb.potcar(
			functional=functional,
		),
		incar = incar,
		kpoints = kpoints,

		# additional files
		metadata = json.dumps(metadata),
	)

	ws.write_input(path)

def make_study_special_files(path, output_dirs):
	# Script for running VASP on all of the trials
	with open(os.path.join(path, SCRIPT_FILENAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

	# File containing list of trial directories, one per line
	with open(os.path.join(path, TRIAL_DIR_LIST_FILENAME), 'w') as f:
		f.write('\n'.join(output_dirs))
		f.write('\n')


def make_scale_study_inputs(path, seps_to_try, **kwargs):

	output_dirs = []

	for trial_num, atomic_sep in enumerate(seps_to_try):
		output_dir = os.path.join(path, 'trial_{}'.format(trial_num))
		make_vasp_input(output_dir,
			atomic_sep=atomic_sep,
			**kwargs
		)

		output_dirs.append(os.path.abspath(output_dir))

	make_study_special_files(path, output_dirs)

def make_kpoints_study_inputs(path, divs_to_try, **kwargs):

	output_dirs = []

	for trial_num, kpointdivs in enumerate(divs_to_try):
		output_dir = os.path.join(path, 'trial_{}'.format(trial_num))
		make_vasp_input(output_dir,
			kpointdivs=kpointdivs,
			**kwargs
		)

		output_dirs.append(os.path.abspath(output_dir))

	make_study_special_files(path, output_dirs)

# arguments shared by all studies
def get_shared_study_args(args):
	return {
		'path':            args.outdir,
		'perturb_dist':    0.2,
		'chain_length':    args.n_atoms,
		'layer_sep':       15.0,
		'forbid_symmetry': False,
		'functional':      'PBE',
		'comment':         'Carbon Chain',
	}

parser = argparse.ArgumentParser()
parser.add_argument('n_atoms', type=int)

subparsers = parser.add_subparsers(help='available studies')

#-------------------------------------------------------------------------
scale_parser = subparsers.add_parser('scale')
scale_parser.add_argument('outdir', type=str)
scale_parser.add_argument('start',  type=float)
scale_parser.add_argument('stop',   type=float)
scale_parser.add_argument('nsteps', type=int)
scale_parser.add_argument('--kpoints', nargs=3, type=int, required=True, metavar=('KX','KY','KZ'))

def handle_scale_mode(args):
	make_scale_study_inputs(
		seps_to_try = np.linspace(args.start, args.stop, args.nsteps),
		kpointdivs  = args.kpoints,
		**get_shared_study_args(args)
	)

scale_parser.set_defaults(func = handle_scale_mode)
#-------------------------------------------------------------------------
kpoints_parser = subparsers.add_parser('kpoints')
kpoints_parser.add_argument('outdir', type=str)
kpoints_parser.add_argument('--scale', type=float, required=True)

def handle_kpoints_mode(args):
	make_kpoints_study_inputs(
		atomic_sep  = args.scale,
		divs_to_try = [[i, 1, 1] for i in range(5,100,5)],
		**get_shared_study_args(args)
	)

kpoints_parser.set_defaults(func = handle_kpoints_mode)
#-------------------------------------------------------------------------

args = parser.parse_args(sys.argv[1:])

# Begin processing for the 'mode' specified by the user.
# (this unusual looking technique is detailed in the argparse docs;
#  'func' is a fake argument supplied by each subparser, which contains a callback)
args.func(args)
