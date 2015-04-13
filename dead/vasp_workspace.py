#!/usr/bin/env python3

#
# keeping this file around for the moment as it wasn't committed before I found
# pymatgen and gave up on this.  ParticleList might be worth bringing back,
# perhaps with just a regular, non-list "add_particle" method.  Just passing data
#  directly into Poscar and Structure, I run the risk of getting things out of order
#  or forgetting to insert something

import math
import sys
import os

import util
import potcar_locator
from workspace_maker import WorkspaceMaker

from enum import Enum, unique

# TODO figure out what this actually is
CARBON_SEP_ANGSTROMS = 1.42

class KPoints:
	pass

# a hashable, immutable type specifying an atomic element
class Species:
	def __init__(self, chemical_symbol):
		self._symbol = chemical_symbol
	def chemical_symbol(self):
		return self._symbol
	def __hash__(self):
		return hash(self._symbol)
	def __str__(self):
		return self._symbol
	def __eq__(self, other):
		return self._symbol == other._symbol

@unique
class Coords(Enum):
	direct    = 1
	cartesian = 2

	def poscar_label(self):
		if self == Coords.direct:
			return 'Direct'
		elif self == Coords.cartesian:
			return 'Cartesian'
		else:
			raise RuntimeError('Logic error: incomplete switch')

@unique
class VelocityMode(Enum):
	manual    = 1
	automatic = 2

class ParticleList:
	def __init__(self, velocity_mode, position_coords=Coords.direct, velocity_coords=Coords.direct):
		expect_type(velocity_mode, VelocityMode)

		self._species_order = []
		self._potcars    = dict()
		self._positions  = dict()
		self._dynamics   = dict()
		self._velocities = dict()

		self._velocity_mode = velocity_mode
		self.set_position_coord_system(position_coords)
		self.set_velocity_coord_system(velocity_coords)

	#----------------

	def set_velocity_mode(self, mode):
		if self._count_particles() > 0:
			raise RuntimeError('Cannot change velocity mode after adding particles')

		expect_type(mode, VelocityMode)
		self._velocity_mode = mode

	def has_velocity(self):
		return self._velocity_mode is VelocityMode.manual

	def _count_particles(self):
		return sum(len(v) for v in self._positions.values())

	#----------------

	def set_position_coord_system(self, coords):
		expect_type(coords, Coords)
		self._position_coords = coords

	def set_velocity_coord_system(self, coords):
		expect_type(coords, Coords)
		self._velocity_coords = coords

	def position_coord_system(self):
		return self._position_coords

	def velocity_coord_system(self):
		return self._velocity_coords

	#----------------

	def add_species(self, species, potcar):
		if species in self._positions:
			raise KeyError('Species {} already added!'.format(repr(species)))

		self._species_order.append(species)

		self._potcars[species]    = potcar
		self._positions[species]  = []
		self._dynamics[species]   = []
		self._velocities[species] = []

	def species(self):
		return iter(self._species_order)

	def species_potcar(self, species):
		return self._potcars[species]

	#----------------

	def add_particles(self, species, position_iter, dynamics_iter, velocity_iter=None):
		if (velocity_iter is None) == self.has_velocity():
			raise ValueError('velocities were provided when none were expected! (or vice versa)')

		ps = list(map(tuple, position_iter))
		ds = list(map(tuple, dynamics_iter))
		vs = list(map(tuple, velocity_iter)) if velocity_iter else None

		if len(ps) != len(ds) or (self.has_velocity() and len(ps) != len(vs)):
			raise ValueError('Non-matching list lengths!')

		self._positions[species].extend(ps)
		self._dynamics[species].extend(ds)

		if self.has_velocity():
			self._velocities[species].extend(vs)

	def species_positions(self, species):
		return list(self._positions[species])

	def species_dynamics(self, species):
		return list(self._dynamics[species])

	def species_velocities(self, species):
		if not self.has_velocity():
			raise RuntimeError("called species_velocities on a ParticleList that doesn't have velocity")
		return list(self._velocities[species])



class LatticeCell:
	def __init__ (self):
		self._raw_scale_factor  = 1.
		self._lattice_vectors   = [
			[1., 0., 0.],
			[0., 1., 0.],
			[0., 0., 1.],
		]

	def set_scale_factor(self, val):
		if not val > 0.0:
			raise ValueError('Scale must be positive (received {})'.format(val))

		self.set_raw_scale_factor(val)

	# for interacting with the raw, vasp-format scale factor
	# (which treats negative numbers as the total cell volume)
	def set_raw_scale_factor(self, val):
		self._raw_scale_factor = val

	def raw_scale_factor(self):
		return self._raw_scale_factor

	def set_lattice_vectors(self, vs):
		if len(vs) != 3:
			raise ValueError("Expected 3 values (received {})".format(len(v)))
		if any(len(v) != 3 for v in vs):
			raise ValueError("Lattice vectors must be of length 3, not {}".format(len(v)))

		self._lattice_vectors = tuple(tuple(v) for v in vs)

	def lattice_vectors(self):
		return self._lattice_vectors


# VASP's INCAR format can only be described as complete and utter mayhem:
#   http://cms.mpi.univie.ac.at/vasp/guide/node91.html
#
# The reader is invited to attempt to reconcile this line from the example:
#   IALGO  =     18    algorithm   NELM   =     60;   NELMIN= 0; NELMDL=  3    # of ELM steps m
# with the description of the format given on the page.
#
# Oh, another fun exercise:  Is the following line a comment, or does it do something?
#   # ISTART = 0; ICHARG = 2
#
# I can do this all day.

class IncarLogical:
	def __init__(self, value):
		self._val = bool(value)
	def value(self):
		return self._val
	def stringify(self):
		return '.TRUE.' if self._val else '.FALSE.'

class IncarReal:
	def __init__(self, value):
		self._val = float(value)
	def value(self):
		return self._val
	def stringify(self):
		# Python includes '.0' for integral floats, and always prints enough precision
		#  to ensure x == float(str(x)). This behavior seems good enough
		return str(self._val)

class IncarInteger:
	def __init__(self, value):
		self._val = int(value)
	def value(self):
		return self._val
	def stringify(self):
		return str(self._val)

class IncarString:
	def __init__(self, value):
		# I'm tempted to reject any string containing a hash, semicolon, equals sign,
		# or a trailing slash almost solely based on how confusing the vasp docs are.
		self._val = str(value)
	def value(self):
		return self._val
	def stringify(self):
		# Deliberately not using repr() because INCAR's string properties don't have quotes.
		return str(self._val)

class Incar:
	def __init__(self, title):
		self._props = {}
		self.set_title(title)

	def set_raw_property_float(self, key, value):
		self._props[key] = IncarReal(value)
	def set_raw_property_bool(self, key, value):
		self._props[key] = IncarLogical(value)
	def set_raw_property_int(self, key, value):
		self._props[key] = IncarInteger(value)
	def set_raw_property_string(self, key, value):
		self._props[key] = IncarString(value)

	def raw_property(self, key):
		return self._props[key].value()

	def raw_property_incar_line(self, key):
		return '{} = {}\n'.format(key, self._props[key].stringify())

	def raw_property_keys(self):
		return self._props.keys()

	def set_title(self, s):
		self.set_raw_property_string('SYSTEM', s)


class VaspWorkspace:
	def __init__(self):
		self._title = 'Untitled'
		self._lattice   = LatticeCell()
		self._particles = ParticleList(VelocityMode.automatic)
		self._kpoints   = KPoints()
		self._incar     = Incar(self._title)

		self._particles_already_added = False

	#-----------
	# methods with special implementation

	def set_title(self, s):
		self._title = s
		self._incar.set_title(s)

	def title(self):
		return self._title

	def species(self):
		return self._particles.species()

	def set_incar_raw(self, key, value):
		if isinstance(value, int):
			self.set_raw_property_int(key, value)
		elif isinstance(value, bool):
			self.set_raw_property_bool(key, value)
		elif isinstance(value, str):
			self.set_raw_property_string(key, value)
		elif isinstance(value, float):
			self.set_raw_property_float(key, value)
		else:
			raise TypeError('Not sure what to do with type {}!'.format(type(value)))

	def set_coord_system(self, position_coords, velocity_coords=None):
		if velocity_coords is None:
			velocity_coords = position_coords
		self._particles.set_position_coord_system(position_coords)
		self._particles.set_velocity_coord_system(velocity_coords)

	def add_particles(self, species, positions, dynamics, velocities=None):
		self._particles.add_particles(species, position_iter=positions, dynamics_iter=dynamics, velocity_iter=velocities)

	#-----------
	# methods delegated to a single member

	# try to keep this list lightweight and manageable

	def add_species(self, *a, **k):   self._particles.add_species(*a, **k)

	def set_scale_factor(self, *a, **k):    self._lattice.set_scale_factor(*a, **k)
	def set_lattice_vectors(self, *a, **k): self._lattice.set_lattice_vectors(*a, **k)

	#----------------

	def make_workspace(self, path):
		wm = WorkspaceMaker()
		wm.bind_writer('POSCAR',  write_poscar_impl,  self._title, self._lattice, self._particles)
		wm.bind_writer('INCAR',   write_incar_impl,   self._incar)
		wm.bind_writer('POTCAR',  write_potcar_impl,  self._particles)
		wm.bind_writer('KPOINTS', write_kpoints_impl, self._kpoints)
		wm.make_workspace(path)


def write_poscar_impl(f, title, lattice, particles):
	expect_type(lattice,   LatticeCell)
	expect_type(particles, ParticleList)

	f.write('{}\n'.format(title))

	f.write(' {:-.15f}\n'.format(lattice.raw_scale_factor()))

	for v in lattice.lattice_vectors():
		f.write('  '.join('{:-.15f}'.format(x) for x in v))
		f.write('\n')

	species_list   = list(particles.species())
	positions_list = list(map(particles.species_positions, species_list))
	dynamics_list  = list(map(particles.species_dynamics,  species_list))

	f.write(' '.join(x.chemical_symbol() for x in species_list))
	f.write('\n')

	f.write(' '.join(str(len(v)) for v in positions_list))
	f.write('\n')

	f.write('Selective dynamics\n')
	f.write('{}\n'.format(particles.position_coord_system().poscar_label())) # Cartesian or Direct

	for (ps, ds) in zip(positions_list, dynamics_list):

		# good point to double check (as zip will simply toss any excess elements)
		assert len(ps) == len(ds)

		for (p, d) in zip(ps, ds):
			pstrs = ['{:-.15f}'.format(x) for x in p]
			dstrs = ['T' if x else 'F' for x in d]
			f.write(' '.join(pstrs + dstrs))
			f.write('\n')

	if particles.has_velocity():
		assert False # FIXME test this crap

		f.write('{}\n'.format(particles.velocity_coord_system().poscar_label())) # Cartesian or Direct

		velocities_list = list(map(particles.species_velocities, species_list))
		for vs in velocities_list:
			for v in vs:
				f.write(' '.join('{:15f}'.format(x) for x in v))
				f.write('\n')

#  TODO: move me
def cubic_lattice(xdim, ydim, zdim):
	return (
		(xdim,  0.0,  0.0),
		( 0.0, ydim,  0.0),
		( 0.0,  0.0, zdim),
	)

def write_kpoints_impl(f, kpoints):
	# FIXME
#	expect_type(kpoints, KPoints)
	f.write("hi i'm kpoints\n")
	f.write(' 0\n') # 0 is automatic generation
	f.write('Gamma point\n') # monkhurst or gamma point
	f.write('9 9 1\n') # mesh parameter (# intersections along each axis)
	f.write('0 0 0\n') # shift


def write_potcar_impl(f, particles):
	expect_type(particles, ParticleList)

	for potcar in map(particles.species_potcar, particles.species()):
		util.validate_file_exists(potcar.path())

		with open(potcar.path()) as fin:
			f.write(fin.read())

def write_incar_impl(f, incar):
	expect_type(incar, Incar)

	for k in incar.raw_property_keys():
		f.write(incar.raw_property_incar_line(k))

def expect_type(value, clazz):
	if not isinstance(value, clazz):
		raise RuntimeError('Expected instance of {} (recieved {})'.format(clazz.__name__, type(value).__name__))

def main():
	make_workspace(sys.argv[1])

if __name__ == '__main__':
	main()
