
from enum import *

import numpy as np
import pymatgen as mg

from .util import validate_type

@unique
class VelocityMode(Enum):
	manual    = 1  # List contains velocities
	automatic = 2  # List does not have velocities (let VASP randomize them)

@unique
class Coords(Enum):
	direct    = 1  # Fractional coords w.r.t. lattice vectors
	cartesian = 2  # Absolute cartesian coords (in whatever units VASP expects)

	def poscar_label(self):
		if self == Coords.direct:
			return 'Direct'
		elif self == Coords.cartesian:
			return 'Cartesian'
		else:
			raise RuntimeError('Logic error: incomplete switch')

# Is intended to help with numerous aspects of building POSCAR and POTCAR:
#  * Keeping species order consistent
#  * Keeping positions, dynamics, and velocities lists in sync
#     (PoxcarBuilder takes all of the data for a particle at once)
#  * Keeping positions together with the lattice vectors that define them,
#     as mixing/matching them only really seems to lead to problems.
#     (pymatgen's Structure already does this, but it lacks the other listed aspects)
class PoxcarBuilder:

	def __init__(self,
		unique_species,
		lattice,
		velocity_mode,
		position_coords = Coords.direct,
		velocity_coords = Coords.direct,
	):

		self._symbols = list(unique_species)
		if len(self._symbols) != len(set(self._symbols)):
			raise ValueError('Species list contains duplicates!')

		validate_type(velocity_mode, VelocityMode)
		self._velocity_mode = velocity_mode

		self._positions  = {x:[] for x in self._symbols}
		self._dynamics   = {x:[] for x in self._symbols}
		self._velocities = {x:[] for x in self._symbols}

		self.set_lattice_vectors(lattice)

		self.set_position_coord_system(position_coords)
		self.set_velocity_coord_system(velocity_coords)

	#----------------

	def set_velocity_mode(self, mode):
		# Protect the "num velocities == num positions" invariant
		if self._count_particles() > 0:
			raise RuntimeError('Cannot change velocity mode after adding particles')

		validate_type(mode, VelocityMode)
		self._velocity_mode = mode

	def has_velocity(self):
		return self._velocity_mode is VelocityMode.manual

	def _count_particles(self):
		return sum(len(v) for v in self._positions.values())

	#----------------

	def set_lattice_vectors(self, lattice):
		arr = np.array(lattice, dtype=float)

		if arr.shape != (3,3):
			raise ValueError('Expected lattice to be an array of shape (3,3) (got {})'.format(self._lattice.shape))

		# detect negative volumes because I'm a bit worried about how pymatgen handles them
		if np.linalg.det(arr) < 0:
			warnings.warn(
				'Volume of lattice cell is negative!  To avoid a bug in VASP, the lattice '
				'vectors may be flipped in the output POSCAR file.  Be aware that this may '
				'INVERT the parity of the structure.'
			) # TODO investigate this (is it really a cause for concern?)

		self._lattice = arr

	#----------------

	def set_position_coord_system(self, coords):
		validate_type(coords, Coords)
		self._position_coords = coords

	def set_velocity_coord_system(self, coords):
		validate_type(coords, Coords)
		self._velocity_coords = coords

	def position_coord_system(self):
		return self._position_coords

	def velocity_coord_system(self):
		return self._velocity_coords

	#----------------

	# For now, seems best to assume all species are specified at initialization.
#	def add_species(self, species, potcar):
#		if species in self._positions:
#			raise KeyError('Species {} already present!'.format(repr(species)))
#
#		self._symbols.append(species)
#
#		self._potcars[species]    = potcar
#		self._positions[species]  = []
#		self._dynamics[species]   = []
#		self._velocities[species] = []

	def species(self):
		return iter(self._symbols)

	#----------------

	def add_particle(self, species, position, dynamics, velocity=None):
		if (velocity is None) == self.has_velocity():
			raise ValueError('velocities were provided when none were expected! (or vice versa)')

		self._positions[species].append(position)
		self._dynamics[species].append(dynamics)

		if self.has_velocity():
			self._velocities[species].append(velocity)

	# Not very useful as it forces you to build dynamics/velocity lists independently of the
	#  positions list, which is the very thing this class was built to prevent.
#	def add_particles(self, species, position_iter, dynamics_iter, velocity_iter=None):
#		if (velocity_iter is None) == self.has_velocity():
#			raise ValueError('velocities were provided when none were expected! (or vice versa)')
#
#		ps = list(map(tuple, position_iter))
#		ds = list(map(tuple, dynamics_iter))
#		vs = list(map(tuple, velocity_iter)) if velocity_iter else None
#
#		if len(ps) != len(ds) or (self.has_velocity() and len(ps) != len(vs)):
#			raise ValueError('Non-matching list lengths!')
#
#		self._positions[species].extend(ps)
#		self._dynamics[species].extend(ds)
#
#		if self.has_velocity():
#			self._velocities[species].extend(vs)

	# TODO: Fix comment
	# TODO: The reason that these methods are the only available accessors is to force the client to iterate
	#        over .species(). (this is the easiest way to guarantee everything is ordered consistently)
	#
	#       I should investigate whether this is still necessary, now that I am aware that VASP allows species
	#        names to be provided at various places in the POSCAR.
#	def species_positions(self, species):
#		return list(self._positions[species])
#
#	def species_dynamics(self, species):
#		return list(self._dynamics[species])
#
#	def species_velocities(self, species):
#		if not self.has_velocity():
#			raise RuntimeError("called species_velocities on a ParticleList that doesn't have velocity")
#		return list(self._velocities[species])

	def _unique_species_list(self):
		return list(self.species())

	# Methods to produce full lists of various properties for all of the atoms as
	#  flattened lists of N elements (ordered by species)

	def _site_species_list(self):
		arr = []
		for species in self.species():
			sublist = [species] * len(self._positions[species])
			arr.extend(sublist)
		return arr

	def _site_position_list(self):
		arr = []
		for species in self.species():
			arr.extend(self._positions[species])
		return arr

	def _site_dynamics_list(self):
		arr = []
		for species in self.species():
			arr.extend(self._dynamics[species])
		return arr

	def _site_velocity_list(self):
		arr = []
		for species in self.species():
			arr.extend(self._velocities[species])
		return arr

	# Output methods

	# TODO: Consider:  Worth forwarding **kwargs to Poscar/Potcar constructor?
	def poscar(self, comment='Automatically generated POSCAR'):
		structure = mg.Structure(
			self._lattice,
			self._site_species_list(),
			self._site_position_list(),
		)

		return mg.io.vaspio.Poscar(
			structure,
			comment,
			selective_dynamics = self._site_dynamics_list(),
			velocities = self._site_velocity_list() if self.has_velocity() else None,
		)

	def potcar(self, functional):
		return mg.io.vaspio.Potcar(
			self._unique_species_list(),
			functional = functional,
		)

