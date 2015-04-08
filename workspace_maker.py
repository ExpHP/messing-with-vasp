#!/usr/bin/env python3
# python 2 sucks and you know it

# 2015 Apr 8

import os

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

		if not os.path.exists(head):
			raise RuntimeError("directory '{}' does not exist".format(head))
		if not os.path.isdir(head):
			raise RuntimeError("'{}' is not a directory".format(head))
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


