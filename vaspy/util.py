
import os

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

def validate_type (value, clazz):
	if not isinstance(value, clazz):
		raise RuntimeError('Expected instance of {} (recieved {})'.format(clazz.__name__, type(value).__name__))
