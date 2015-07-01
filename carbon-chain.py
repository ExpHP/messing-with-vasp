
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

def main(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('n_atoms', type=int)

	add_subparsers(parser)
	args = parser.parse_args(sys.argv[1:])

	if hasattr(args, 'func'): # workaround for Python 3.3+: http://bugs.python.org/issue16308
		# Begin processing for the 'mode' specified by the user.
		# (this unusual looking technique is detailed in the argparse docs;
		#  'func' is a fake argument supplied by each subparser, which contains a callback)
		args.func(args)
	else:
		parser.print_usage(file=sys.stderr)
		print('Error: No mode specified', file=sys.stderr)
		sys.exit(1)

def add_subparsers(parser):
	subparsers = parser.add_subparsers(help='available studies')

	add_scale_subparser(subparsers, 'scale')
	add_kpoints_subparser(subparsers, 'kpoints')

def add_scale_subparser(subparsers, command):
	sub = subparsers.add_parser(command)
	sub.add_argument('outdir', type=str)
	sub.add_argument('start',  type=float)
	sub.add_argument('stop',   type=float)
	sub.add_argument('nsteps', type=int)
	sub.add_argument('--kpoints', type=int, required=True, metavar='KZ')

	def callback(args):
		make_scale_study_inputs(
			seps_to_try = np.linspace(args.start, args.stop, args.nsteps),
			kpointdivs  = [1, 1, args.kpoints],
			**get_shared_study_args(args)
		)
	sub.set_defaults(func = callback)

def add_kpoints_subparser(subparsers, command):
	sub = subparsers.add_parser(command)
	sub.add_argument('outdir', type=str)
	sub.add_argument('--scale', type=float, required=True)

	def callback(args):
		make_kpoints_study_inputs(
			atomic_sep  = args.scale,
			divs_to_try = [[1, 1, x] for x in range(5,100,5)],
			**get_shared_study_args(args)
		)
	sub.set_defaults(func = callback)

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

#-------------------------------------------------------------------------

def carbon_chain_poxcar (
	chain_length,    # number of links to include of the chain
	atomic_sep,      #
	layer_sep,       # separation between periodic copies on both axes
	specie = 'C',
):
	dimensions = np.array([
		layer_sep,
		layer_sep,
		chain_length * atomic_sep,
	], dtype=float)

	lattice = np.eye(3) * dimensions

	pb = PoxcarBuilder([specie], lattice, VelocityMode.automatic)

	for x in periodic_points(0., 1., chain_length):
		pb.add_particle(specie, [1/2, 1/2, x], [True, True, True])

	return pb

#-------------------------------------------------------------------------

# makes a study where an argument to make_vasp_input is varied
def make_study_inputs(path, argname, values, **kwargs):

	output_dirs = []

	# Trials
	for trial_num, value in enumerate(values):
		output_dir = os.path.join(path, 'trial_{}'.format(trial_num))
		output_dirs.append(os.path.abspath(output_dir))

		kwargs[argname] = value # vary the argument
		make_vasp_input(output_dir, **kwargs)

	# Script for running VASP on all of the trials
	with open(os.path.join(path, SCRIPT_FILENAME),'w') as f:
		f.write(SCRIPT_CONTENTS)

	# File containing list of trial directories
	with open(os.path.join(path, TRIAL_DIR_LIST_FILENAME), 'w') as f:
		f.write('\n'.join(output_dirs))
		f.write('\n')

def make_scale_study_inputs(path, seps_to_try, **kwargs):
	make_study_inputs(path, 'atomic_sep', seps_to_try, **kwargs)

def make_kpoints_study_inputs(path, divs_to_try, **kwargs):
	make_study_inputs(path, 'kpointdivs', divs_to_try, **kwargs)

def make_vasp_input(path, *, kpointdivs, comment, forbid_symmetry, functional, perturb_dist, **kwargs):
	assert hasattr(kpointdivs, '__iter__')
	assert isinstance(comment, str)
	assert isinstance(forbid_symmetry, bool)
	assert isinstance(perturb_dist, float)

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

	metadata = {
		'kpointdivs':kpointdivs,
		'functional':functional,
		'perturb_dist':perturb_dist,
		'forbid_symmetry':forbid_symmetry,
	}

	# FIXME FIXME FIXME
	# Here, we insert all structure parameters into the metadata.  The purpose is because
	#  oftentimes those parameters are the ones we're most interested in (such as scale).
	#
	# However, this is all kinds of bad:
	#  * Not all structure arguments are necessarily desirable to have in the metadata.
	#  * What if we want a structure argument to be of a non JSON-encodable type?
	#  * It creates a dependency between names of function arguments in code, and the output file.
	# In other words, madness.
	metadata = dict_union(metadata, kwargs)

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

#-------------------------------------------------------------------------

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


# Fail-fast alternative to dict.update.
# Combines two dicts, throwing an error if any keys are shared.
def dict_union(d1, d2):
	shared_keys = set(d1) & set(d2)
	if len(shared_keys) != 0:
		key = shared_keys.pop()
		msg = 'duplicate key {} (values: {} vs {})'.format(repr(key), repr(d1[key]), repr(d2[key]))
		msg += '' if len(shared_keys) == 0 else ' (and {} more)'.format(len(shared_keys))
		raise ValueError(msg)

	result = dict(d1)
	result.update(d2)
	return result

#-------------------------------------------------------------------------

if __name__ == '__main__':
	main(sys.argv)

