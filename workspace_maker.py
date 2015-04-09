#!/usr/bin/env python3
# python 2 sucks and you know it

# 2015 Apr 8

import os

# TODO: Wrap the use of environment vars in a more future-friendly config-providing class or module

def validate_file_exists (path, msg=''):
	if not os.path.exists(path):
		raise RuntimeError('"{}" does not exist. '.format(path) + msg)
	if not os.path.isfile(path):
		raise RuntimeError('"{}" is not a file. '.format(path) + msg)

def validate_dir_exists (path, msg=''):
	if not os.path.exists(path):
		raise RuntimeError('"{}" does not exist. '.format(path) + msg)
	if not os.path.isdir(path):
		raise RuntimeError('"{}" is not a directory. '.format(path) + msg)

# Abstract base class for creating a VASP workspace, to eliminate any boring
#  filesystem-related copy pasta.
class VaspWorkspaceMaker:

	# Public base class API

	def make_workspace (self, relpath):
		workspace_path = os.path.abspath(relpath)
		self._validate_destination(workspace_path)
		self._create_workspace(workspace_path)

	# Abstract methods - IMPLEMENT THESE ON BASE CLASS

	def write_poscar (self, f):
		raise NotImplementedError("write_poscar not implemented on " + type(self).__name__)
	def write_incar (self, f):
		raise NotImplementedError("write_incar not implemented on " + type(self).__name__)
	def write_potcar (self, f):
		raise NotImplementedError("write_potcar not implemented on " + type(self).__name__)
	def write_kpoints (self, f):
		raise NotImplementedError("write_kpoints not implemented on " + type(self).__name__)

	# Helper methods

	# NOTE: Not actually sure if this logic really belongs here since I might want e.g. a --force
	#       command line option, and it wouldn't make sense for this class to handle that... :/
	def _validate_destination (self, workspace_path):
		head,tail = os.path.split(workspace_path)

		validate_dir_exists(head)

		if os.path.exists(workspace_path):
			raise RuntimeError("'{}' already exists".format(workspace_path))

	def _create_workspace (self, path):
		os.mkdir(path)

		with self._open_output_file(path, 'POSCAR') as f:
			self.write_poscar(f)

		with self._open_output_file(path, 'INCAR') as f:
			self.write_incar(f)

		with self._open_output_file(path, 'POTCAR') as f:
			self.write_potcar(f)

		with self._open_output_file(path, 'KPOINTS') as f:
			self.write_kpoints(f)

	def _open_output_file (self, workspace, fname):
		return open(os.path.join(workspace, fname), 'w')


#--------------------------------------------------------

#--------------------------------------------------------
# TODO: This static code really shouldn't be in this file as its presence makes the
#       abstract base class implicitly dependent on the implementation
#       (the ABC cannot be imported if the locpot path doesn't exist)
ENV_POTCAR_DIR = 'VASP_POTCAR_DIR'

MSG_BAD_POTCAR_DIR = "Please link this path (or set the environment variable {}) to point to Vasp's potentials/potcar directory".format(ENV_POTCAR_DIR)

POTCAR_DIR = os.environ.get(ENV_POTCAR_DIR)
if POTCAR_DIR is None:
	POTCAR_DIR = os.path.join(os.getcwd(), 'locpot') # default setting

validate_dir_exists(POTCAR_DIR, MSG_BAD_POTCAR_DIR)

#--------------------------------------------------------
# Helper methods which partially implement writing the files

def write_composite_potcar(f, names_iter):
	for name in names_iter:
		write_vasp_builtin_potcar(f, name)

# Writes the contents of a vasp builtin potcar file into the provided
#  output file handle
def write_vasp_builtin_potcar(fout, name):
	src_potcar_path = os.path.join(POTCAR_DIR, name, 'POTCAR')

	validate_file_exists(src_potcar_path)
	with open(src_potcar_path) as fin:
		s = fin.read()

	fout.write(s)

